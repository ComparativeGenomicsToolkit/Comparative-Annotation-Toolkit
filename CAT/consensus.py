"""
Generates consensus gene set.

This module takes as input the genePreds produced by transMap, AugustusTM(R) and AugustusCGP and generates a consensus
of these, producing a filtered gene set.

This process relies on a combination of metrics and evaluations loaded to a sqlite database by the classify module.

Transcript scoring functions:

structure score = 100 *  (0.5 * (1 - # missing introns / # parent introns)
                          + 0.45 * (1 - # missing exons / # parent exons)
                          + 0.05 * # exon gain / # parent exons)
This function is a weighted average of the structural changes seen in the alignments. Missing original introns is
weighted slightly higher than losing exons which is slightly higher than gaining exons because gaining exons may be
a real change. We multiply by 100 to put this score on the same scale as the base transcript score.
# parent introns/exons will be adjusted for CDS alignments to only include coding introns/exons.
Exists on the range(0, 100)

evaluation penalties = 2 * I(in frame stop)
                         + I(coding indel)
                         + I(CdsStartStat = imcpl and == cmpl in reference)
                         + I(CdsEndStat = imcpl and == cmpl in reference)
This function uses classifiers in a binary fashion as indicator functions to penalize the final transcript score for
problems. Only applies to coding transcripts. This function does not penalize transcripts who are messed up in the
reference, as this would bias towards incorrect Augustus corrections.
Exists in the range (0, 5) where 0 means no penalties and 5 is the maximum penalty.

evaluation bonus = 2 * I(CdsStartStat == cmpl) + 2 * I(CdsEndStat == cmpl)
This function looks at the CdsStartStat/CdsEndStat and provides up to 2 bonus points for complete CDS ends.
Exists in [0, 2, 4]


consensus score = 50 * (1 - badness) + 25 * structure score + 25 * # supported junctions
Exists in the range (0, 100). Higher is better.

coding consensus score = consensus score + cds bonus + evaluation bonus - evaluation penalties
For coding transcripts consensus.py modifies the consensus score by removing as much as 5 and adding as much as 4

Choosing CDS vs mRNA:

align_transcripts.py aligns transMap/AugustusTM(R) transcripts both in both CDS space and mRNA space.
For CGP transcripts, only CDS space alignments are performed. Thus, we do not have a fair
comparison if we use the mRNA alignments. Since CDS sequences are more conserved, it can be expected that if alignments
are evaluated all together that the CDS alignments will have the highest scores unless something has gone seriously
wrong (often, incomplete mapping particularly for genes with long UTRs). In those cases, the mRNA alignment will win.
For this reason, all scores are considered simultaneously.

The filter requirements below ensure that the results are at least of sufficient quality.

Filtering transcripts:

Before scoring, the transcripts are filtered based on a series of minimum cutoffs. These are:
AlnIdentity: mRNA: >70% CDS: >80%
AlnCoverage: mRNA: >50% CDS: >90%
Intron Inequality: 2 * # missing original introns <= # reference introns - 1 or # reference introns < 5
PercentUnknownBases: mRNA: <5% CDS: <2%

These filters throw out low quality alignments via the identity/coverage filter as well as the percent unknown bases.
The intron inequality is an important filter to preventing retroposed pseudogenes from being assigned as the parent.

If the filter process throws out all transcripts for a gene, then one transcript with the lowest badness score
will be kept to represent the locus.

Paralogy: If a transMap transcript mapped to more than one place, then we make use of the synteny score to decide
which locus is more correct. Generally, paralogous alignments are one of two things: actual paralogy, and alignment
chain breaking rearrangements that lead to fractured transcripts. In the first case, the synteny score should help. In
the second case, the synteny score will be generally the same. If this happens, then hopefully AUGUSTUS in one of its
forms rescued the transcript.

"""
import collections
import pandas as pd
import tools.mathOps
import tools.transcripts
import tools.sqlInterface
import tools.intervals

cds_cutoffs = {'AlnIdentity': 0.8, 'AlnCoverage': 0.9, 'PercentUnknownBases': 0.02}
mrna_cutoffs = {'AlnIdentity': 0.7, 'AlnCoverage': 0.5, 'PercentUnknownBases': 0.05}


def consensus(args):
    """
    Entry point for consensus finding module.
    :param args: Argument namespace from luigi
    """
    # load all genePreds
    tx_dict = tools.transcripts.load_gps(args.gp_list)
    ref_tx_dict = tools.transcripts.get_gene_pred_dict(args.ref_gp)

    # load database tables
    ref_df = tools.sqlInterface.load_reference(args.ref_genome_db)
    aln_eval_df = tools.sqlInterface.load_alignment_evaluation(args.db_path)

    # load alignment-mode specific tables, scoring the alignments
    aln_mode_dfs = {}
    for aln_mode in ['CDS', 'mRNA']:
        intron_df = tools.sqlInterface.load_intron_vector(args.db_path, aln_mode, args.transcript_modes, tx_dict)
        tgt_df = tools.sqlInterface.load_classifications(args.db_path, aln_mode, args.transcript_modes, ref_tx_dict)
        merged_df = merge_ref_tgt(ref_df, tgt_df, intron_df)
        filtered_df = filter_alignments(merged_df, aln_mode)
        scored_df = score_df(filtered_df)
        aln_mode_dfs[aln_mode] = scored_df

    # metrics will hold values during consensus finding that will be plotted afterwards
    metrics = {}

    # resolve potential conflicts before consensus finding
    alns_to_remove = set()
    if args.resolve_split_genes is True:
        alns_to_remove.update(resolve_split_genes(aln_eval_df, ref_df, tx_dict))
        metrics['Genes Lost To Split Chromosomes'] = num_removed

    # remove filtered alignments



    # begin consensus finding


###
# Filtering functions
###


def filter_alignments(merged_df, aln_mode):
    """
    Takes the combined DataFrame produced by load_classifications() and filters it based on the filter criteria
    Applies a filter for coverage, identity, unknown base count, and the intron inequality
    :param merged_df: DataFrame produced by load_classifications()
    :param aln_mode: one of ('CDS', 'mRNA')
    :return: DataFrame
    """
    cutoffs = cds_cutoffs if aln_mode == 'CDS' else mrna_cutoffs
    return merged_df[(merged_df['AlnIdentity'] > cutoffs['AlnIdentity']) &
                     (merged_df['AlnCoverage'] > cutoffs['AlnCoverage']) &
                     (merged_df['PercentUnknownBases'] <= cutoffs['PercentUnknownBases']) &
                     ((2 * merged_df['NumMissingIntrons'] <= merged_df['NumReferenceIntrons'] - 1)
                      | (merged_df['NumReferenceIntrons'] < 5))]


###
# Merging functions
###


def merge_ref_tgt(ref_df, tgt_df, intron_df):
    """
    Merges the two dataframes into one using what is effectively a join statement. Reorganizes the index such that
    we now have a hierarchy of gene -> transcripts -> alignments
    For columns which share the same name (StartCodon/StopCodon), we append a unique identifier
    :param ref_df: DataFrame from sqlInterface.load_reference()
    :param tgt_df: DataFrame from sqlInterface.load_classifications()
    :return: DataFrame
    """
    df = pd.merge(tgt_df.reset_index(), ref_df.reset_index(), on='TranscriptId', how='inner', suffixes=['_Tgt', '_Ref'])
    df = pd.merge(df, intron_df.reset_index(), on='AlignmentId', how='inner')
    return df.set_index(['GeneId', 'TranscriptId', 'AlignmentId'])


###
# Scoring functions
###


def score_df(df):
    """
    Scores a merged and filtered DataFrame produced by filter_alignments()
    For now, computing structure score a second time and storing. This is for debug purposes mostly.
    :param df: DataFrame
    :return: DataFrame
    """
    scores = collections.defaultdict(list)
    for _, s in df.iterrows():
        scores['ConsensusScore'].append(consensus_score(s))
        scores['StructureScore'].append(structure_score(s))
    return df.assign(ConsensusScore=scores['ConsensusScore'], StructureScore=scores['StructureScore'])


def consensus_score(s):
    """
    consensus score = 50 * (1 - badness) + 25 * structure score + 25 * # supported junctions
    Exists in the range (0, 100). Higher is better.

    coding consensus score = consensus score + cds bonus + evaluation bonus - evaluation penalties
    For coding transcripts consensus.py modifies the consensus score by removing as much as 5 and adding as much as 4

    Exists in the range (-5, 104)

    :param s: Pandas Series
    :return: float
    """
    score = 50 * (1 - s.Badness) + 25 * structure_score(s) + 25 * supported_junctions(s)
    if s.TranscriptBiotype == 'protein_coding':
        score = score + evaluation_bonus(s) - evaluation_penalties(s)
        assert -5 <= score <= 104
    else:
        assert 0 <= score <= 100
    return score


def indicator(item):
    """Pass a Series member and returns a 1 if the value in the member is > 0"""
    return 1 if item > 0 else 0


def structure_score(s):
    """
    structure score = 0.5 * (1 - # missing introns / # parent introns)
                        + 0.45 * (1 - # missing exons / # parent exons)
                        + 0.05 * # exon gain / # parent exons
    This function is a weighted average of the structural changes seen in the alignments. Missing original introns is
    weighted slightly higher than losing exons which is slightly higher than gaining exons because gaining exons may be
    a real change.
    Parent introns/exons will be adjusted for CDS alignments to only include coding introns/exons by the sqlInterface
    query.
    Exists on the range(0, 1)

    :param s: Pandas Series
    :return: float between 0 and 1
    """
    present_introns = 1 - tools.mathOps.format_ratio(s.NumMissingIntrons, s.NumReferenceIntrons, resolve_nan=1)
    present_exons = 1 - tools.mathOps.format_ratio(s.NumMissingExons, s.NumReferenceExons)
    exon_gain = tools.mathOps.format_ratio(s.ExonGain, s.NumReferenceExons)
    score = 0.5 * present_introns + 0.45 * present_exons + 0.05 * exon_gain
    assert 0 <= score <= 1, s
    return score


def supported_junctions(s):
    """
    Returns the percent of intron junctions supported by RNA-seq as calculated by homGeneMapping
    :param s: Pandas Series
    :return: float between 0 and 1
    """
    percent_supported = tools.mathOps.format_ratio(s.NumSupportedIntrons, s.NumIntrons, resolve_nan=1)
    assert 0 <= percent_supported <= 1, s
    return percent_supported


def evaluation_penalties(s):
    """
    evaluation penalties = 2 * I(in frame stop)
                         + I(coding indel)
                         + I(CdsStartStat = imcpl and == cmpl in reference)
                         + I(CdsEndStat = imcpl and == cmpl in reference)
    This function uses classifiers in a binary fashion as indicator functions to penalize the final transcript score for
    problems. Only applies to coding transcripts. This function does not penalize transcripts who are messed up in the
    reference, as this would bias towards incorrect Augustus corrections.
    Exists in the range (0, 5) where 0 means no penalties and 5 is the maximum penalty.

    These only apply to protein coding transcripts

    :param s: Pandas Series
    :return: integer in range(0, 5)
    """
    # we have to use .get because if somehow we evaluated zero of a evaluation it will not exist in the table
    in_frame_stop_penalty = 2 * indicator(s.get('InFrameStop', 0))
    indel_penalty = indicator(s.get('CodingDeletion', 0) + s.get('CodingInsertion', 0))
    # codon classifiers have 1 for cmpl and 0 otherwise
    start_penalty = 1 if s.StartCodon_Tgt == 0 and s.StartCodon_Ref == 1 else 0
    stop_penalty = 1 if s.StopCodon_Tgt == 0 and s.StopCodon_Ref == 1 else 0
    penalty = in_frame_stop_penalty + indel_penalty + start_penalty + stop_penalty
    assert 0 <= penalty <= 5, s
    return penalty


def evaluation_bonus(s):
    """
    evaluation bonus = 2 * I(CdsStartStat == cmpl) + 2 * I(CdsEndStat == cmpl)
    This function looks at the CdsStartStat/CdsEndStat and provides up to 2 bonus points for complete CDS ends.
    Exists in [0, 2, 4]

    These only apply to protein coding transcripts

    :param s: Pandas Series
    :return: integer in set [0, 2, 4]
    """
    start_bonus = 2 if s.StartCodon_Tgt == 1 else 0
    stop_bonus = 2 if s.StopCodon_Tgt == 1 else 0
    return start_bonus + stop_bonus


###
# Consensus finding functions
###






###
# Helper functions
###

