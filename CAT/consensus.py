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

import tools.intervals
import tools.mathOps
import tools.sqlInterface
import tools.transcripts
import tools.nameConversions

cds_cutoffs = {'AlnIdentity': 0.8, 'AlnCoverage': 0.9, 'PercentUnknownBases': 0.02}
mrna_cutoffs = {'AlnIdentity': 0.7, 'AlnCoverage': 0.5, 'PercentUnknownBases': 0.05}


def consensus(args, consensus_gp):
    """
    Entry point for consensus finding module.
    :param args: Argument namespace from luigi
    """
    # load all genePreds
    tx_dict = tools.transcripts.load_gps(args.gp_list)
    ref_tx_dict = tools.transcripts.get_gene_pred_dict(args.annotation_gp)

    # load reference data
    ref_df = tools.sqlInterface.load_reference(args.ref_db_path)

    # load alignment-mode specific tables, scoring the alignments
    aln_mode_dfs = {}
    for aln_mode in ['CDS', 'mRNA']:
        intron_df = tools.sqlInterface.load_intron_vector(args.db_path, aln_mode, args.transcript_modes, tx_dict)
        tgt_df = tools.sqlInterface.load_classifications(args.db_path, aln_mode, args.transcript_modes, ref_tx_dict)
        merged_df = merge_ref_tgt(ref_df, tgt_df, intron_df)
        filtered_df = filter_alignments(merged_df, aln_mode)
        scored_df = score_df(filtered_df)
        indexed_scored_df = scored_df.set_index(['GeneId', 'TranscriptId'])
        aln_mode_dfs[aln_mode] = indexed_scored_df.sort_index()
        # some of these dataframes need to be saved for later
        if aln_mode == 'mRNA':
            # this will be needed for gene rescue
            unfiltered_df = merged_df.set_index('GeneId')
        if aln_mode == 'CDS':
            cgp_intron_df = intron_df[intron_df.AlignmentId.str.contains('jg')].set_index('AlignmentId')

    # construct a map to iterate over
    gene_transcript_map = tools.sqlInterface.get_gene_transcript_map(args.ref_db_path)

    # metrics will hold values during consensus finding that will be plotted afterwards
    metrics = {}

    # begin consensus finding
    consensus_tx_dict = {}  # stores mapping between target ID and aln_id

    # incorporate every CGP that was not assigned a parental gene
    novel_transcripts = {tx.name: tx.name for tx in tx_dict.itervalues() if 'jg' in tx.name2}
    novel_genes = {tx.name2 for tx in tx_dict.itervalues() if 'jg' in tx.name2}
    metrics['CGP'] = {'Novel genes': len(novel_genes), 'Novel transcripts': len(novel_transcripts)}
    consensus_tx_dict.update(novel_transcripts)

    # main consensus finding loop
    aln_modes = collections.Counter()
    tx_modes = collections.Counter()
    failed_consensus = {'Gene': set(), 'Transcript': set()}
    for gene, tx_list in gene_transcript_map.iteritems():
        gene_seen = False
        for tx in tx_list:
            best = find_best_score(aln_mode_dfs, gene, tx)
            if best is None:
                failed_consensus['Transcript'].add(tx)
            else:
                gene_seen = True
                best_id, tx_mode, aln_mode = best
                consensus_tx_dict[tx] = best_id
                tx_modes[tx_mode] += 1
                aln_modes[aln_mode] += 1
        if gene_seen is False:
            failed_consensus['Gene'].add(gene)

    # save metrics
    metrics['Initial Consensus'] = {'Alignment Modes': aln_modes, 'Transcript Modes': tx_modes,
                                    'Gene Failed': len(failed_consensus['Gene']),
                                    'Transcript Failed': len(failed_consensus['Transcript'])}

    # rescue genes that were lost
    metrics['Gene Rescue'] = collections.Counter()
    for gene in failed_consensus['Gene']:
        # use the unfiltered_df which is the mRNA alignment data without filtering
        rescue = rescue_missing_gene(unfiltered_df, gene)
        if rescue is not None:
            rescued_tx, aln_id = rescue
            consensus_tx_dict[rescued_tx] = aln_id
            metrics['Gene Rescue']['Rescued'] += 1
        else:
            metrics['Gene Rescue']['Failed rescue'] += 1

    # incorporate every CGP as a novel isoform if it has novel supported splice junctions
    novel_isoforms, transcripts_to_include = find_novel_cgp_splices(tx_dict, consensus_tx_dict, cgp_intron_df,
                                                                    gene_transcript_map)
    metrics['CGP']['Novel isoforms'] = novel_isoforms
    consensus_tx_dict.update(transcripts_to_include)

    # deduplicate consensus
    deduplicated_consensus, dup_count = deduplicate_consensus(consensus_tx_dict, tx_dict, unfiltered_df)
    metrics['Duplicates removed'] = dup_count

    # sanity check
    assert len(set(deduplicated_consensus.values())) == len(deduplicated_consensus)

    # write out the results
    with consensus_gp.open('w') as outf:
        for tx_id, aln_id in consensus_tx_dict.iteritems():
            tx = tx_dict[aln_id]
            tx.name = tx_id
            outf.write('\t'.join(tx.get_gene_pred()) + '\n')
    return metrics


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
                     ((2 * merged_df['NumMissingIntrons'] <= merged_df['NumReferenceIntrons'] - 1) |
                      (merged_df['NumReferenceIntrons'] < 5))]


def deduplicate_consensus(consensus_tx_dict, tx_dict, unfiltered_df):
    """
    In the process of consensus building, we may find that we have ended up with more than one transcript for a gene
    that are actually identical. Remove these, picking the best based on their badness, favoring the coding transcript
    """
    def resolve_duplicate(gp_list, df):
        """resolve first by biotype then by score"""
        tx_biotypes = {aln_id: df.ix[aln_id].TranscriptBiotype for tx_id, aln_id in gp_list}
        if len(tx_biotypes) > 1 and 'protein_coding' in tx_biotypes:
            gp_list = [aln_id for aln_id, biotype in tx_biotypes.iteritems() if biotype == 'protein_coding']
            if len(gp_list) == 1:
                return gp_list[0]
        badness_scores = {(tx_id, aln_id): df.ix[aln_id].Badness for tx_id, aln_id in gp_list}
        badness_scores = sorted(badness_scores.iteritems(), key=lambda ((tx_id, aln_id), badness): badness)
        return badness_scores[-1][0]

    # re-index the dataframe
    df = unfiltered_df.reset_index().set_index(['AlignmentId'])

    # build a dictionary mapping duplicates making use of hashing intervals
    duplicates = collections.defaultdict(list)
    for tx_id, aln_id in consensus_tx_dict.iteritems():
        tx = tx_dict[aln_id]
        duplicates[frozenset(tx.exon_intervals)].append([tx_id, aln_id])

    # begin iterating
    deduplicated_consensus = {}
    dup_count = 0
    for gp_list in duplicates.itervalues():
        if len(gp_list) > 1:
            dup_count += 1
            # we have duplicates to collapse - which has the lowest badness?
            tx_id, aln_id = resolve_duplicate(gp_list, df)
            deduplicated_consensus[tx_id] = aln_id
        else:
            tx_id, aln_id = gp_list[0]
            deduplicated_consensus[tx_id] = aln_id
    return deduplicated_consensus, dup_count


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
    df = pd.merge(tgt_df, ref_df, on='TranscriptId', how='inner', suffixes=['_Tgt', '_Ref'])
    return pd.merge(df, intron_df, on='AlignmentId', how='inner')


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


def find_best_score(aln_mode_dfs, gene, tx):
    """
    Extracts the best score for all alignments for a given transcript, handling missing data
    :param aln_mode_dfs: dictionary of DataFrames produced by consensus()
    :param gene: gene ID
    :param tx: tx ID
    :return: best_id, tx_mode, aln_mode or None
    """
    def handle_missing(df, gene, tx):
        try:
            return df.ix[gene, tx]
        except KeyError:
            return None

    def find_best(cds, mrna):
        if cds is None:
            mrna_r = mrna.sort_values('ConsensusScore').iloc[-1]
            return mrna_r.AlignmentId, tools.nameConversions.alignment_type(mrna_r.AlignmentId), 'mRNA'
        elif mrna is None:
            cds_r = cds.sort_values('ConsensusScore').iloc[-1]
            return cds_r.AlignmentId, tools.nameConversions.alignment_type(cds_r.AlignmentId), 'CDS'
        mrna_r = mrna.sort_values('ConsensusScore').iloc[-1]
        cds_r = cds.sort_values('ConsensusScore').iloc[-1]
        if mrna_r.ConsensusScore >= cds_r.ConsensusScore:
            return mrna_r.AlignmentId, tools.nameConversions.alignment_type(mrna_r.AlignmentId), 'mRNA'
        else:
            return cds_r.AlignmentId, tools.nameConversions.alignment_type(cds_r.AlignmentId), 'CDS'

    cds = handle_missing(aln_mode_dfs['CDS'], gene, tx)
    mrna = handle_missing(aln_mode_dfs['mRNA'], gene, tx)
    if cds is None and mrna is None:
        return None
    else:
        return find_best(cds, mrna)


def rescue_missing_gene(unfiltered_df, gene):
    """
    Similar to find_best_score, but resolves on the gene level picking based on badness
    :param unfiltered_df: mRNA alignment mode DataFrame
    :param gene: gene ID to rescue
    :return: best_id or None
    """
    try:
        df = unfiltered_df.ix[gene]
    except KeyError:
        return None

    if not isinstance(df, pd.core.series.Series):
        # if exactly one alignment is being evaluated, it will be a series
        df = df.sort_values('Badness')
        gene_biotype = df.GeneBiotype[0]
        if gene_biotype == 'protein_coding':
            # only consider coding transcripts
            coding_df = df[df.TranscriptBiotype == 'protein_coding']
            if len(coding_df) != 0:
                # edge case where this coding gene has no coding transcripts
                return coding_df.iloc[0].TranscriptId, coding_df.iloc[0].AlignmentId
        return df.iloc[0].TranscriptId, df.iloc[0].AlignmentId
    else:
        return df.TranscriptId, df.AlignmentId


def find_novel_cgp_splices(tx_dict, consensus_tx_dict, cgp_intron_df, gene_transcript_map):
    """
    Finds novel splice junctions in CGP transcripts. If there are any, these get included as a novel isoform.
    :param tx_dict: Combined dictionary of GenePredTranscript objects
    :param consensus_tx_dict: Dict mapping {tx_id: aln_id} for transcripts that are in the consensus set
    :param cgp_intron_df: DataFrame of intron vectors
    :param gene_transcript_map: Dict mapping {gene_id: tx_id}
    :return:
    """
    def create_cgp_dict(consensus_tx_set, tx_dict):
        """create a dict mapping ensembl gene names to associated comparative transcripts"""
        cgp_dict = collections.defaultdict(list)
        for aln_id, tx in tx_dict.iteritems():
            if 'jg' in aln_id and 'jg' not in tx.name2 and aln_id not in consensus_tx_set:
                cgp_dict[tx.name2].append(tx)
        return cgp_dict

    def find_existing_splices(tx_dict, consensus_tx_dict, tx_ids):
        existing_splices = set()
        for tx in tx_ids:
            try:  # this transcript may not be in the consensus
                consensus_tx = tx_dict[consensus_tx_dict[tx]]
            except KeyError:
                continue
            existing_splices.update(consensus_tx.intron_intervals)
        return existing_splices

    def find_supported_intervals(tx, cgp_intron_df):
        intervals = set()
        intron_vector = cgp_intron_df.ix[tx.name].IntronVector
        for interval, intron_score in zip(*[tx.intron_intervals, intron_vector]):
            if intron_score > 0:
                intervals.add(interval)
        return intervals

    novel_isoforms = 0
    transcripts_to_include = {}
    consensus_tx_set = set(consensus_tx_dict.itervalues())
    cgp_dict = create_cgp_dict(consensus_tx_set, tx_dict)
    for gene_id, tx_ids in gene_transcript_map.iteritems():
        existing_splices = find_existing_splices(tx_dict, consensus_tx_dict, tx_ids)
        for tx in cgp_dict[gene_id]:
            tx_supported_intervals = find_supported_intervals(tx, cgp_intron_df)
            if len(tx_supported_intervals - existing_splices) > 0:
                novel_isoforms += 1
                transcripts_to_include[tx.name] = tx.name
    return novel_isoforms, transcripts_to_include
