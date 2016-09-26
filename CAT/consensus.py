"""
Generates consensus gene set.

This module takes as input the genePreds produced by transMap, AugustusTM(R) and AugustusCGP and generates a consensus
of these, producing a filtered gene set.

This process relies on a combination of metrics and evaluations loaded to a sqlite database by the classify module.

Transcript scoring functions:

structure score = 0.5 * (1 - # missing introns / # parent introns) + 0.5 * (1 - # missing exons / # parent exons)
This function is a weighted average of the structural changes seen in the alignments.
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

GFF3 tags generated in this process:
1. source_transcript: The name of the parent transcript, if it exists
2. source_gene: The name of the parent gene, if it exists
3. transcript_mode: The name of the mode of operation that generated this transcript
4. score: the consensus score of the transcript
5. collapsed_duplicate: A comma separated list of alternate IDs for this transcript
6. failed_gene: This transcript is the single representative for a failed transcript
7. novel_sequence:

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


def consensus(args):
    """
    Entry point for consensus finding algorithm. Main consensus finding logic is here.

    Consensus process:
        Loop through genes, then through each transcript. For each transcript, look first at mRNA space alignments.
        Pick the highest. If none pass filtering, look at the CDS space alignments, which may incorporate CGP.
        If no transcripts for a gene pass filtering, resolve based on the lowest CDS badness score, and flag as poor.

    After this process, if CGP was used, the CGP transcripts are evaluated for providing either A) novel genes, or
    B) novel splice junctions.

    :param args: Argument namespace from luigi
    """
    def evaluate_seen(gene_seen, best_row):
        """a gene is only seen if the transcript biotype matches the gene biotype"""
        if gene_seen is False:
            gene_seen = True if best_row.TranscriptBiotype == best_row.GeneBiotype else False
        return gene_seen

    def evaluate_ties(best_rows):
        """Find out how many transcript modes agreed on this"""
        return ','.join(sorted(set([tools.nameConversions.alignment_type(x) for x in best_rows.index])))

    def evaluate_tx_alns(filtered_df, tx_id, aln_mode, gene_seen):
        """evaluates a row of either type of alignment"""
        best_rows = find_best_score(filtered_df, tx_id)
        if best_rows is not None:
            metrics['Transcript Modes'][evaluate_ties(best_rows)] += 1
            aln_id, best_row = best_rows.iterrows().next()
            gene_seen = evaluate_seen(gene_seen, best_row)
            consensus[aln_id] = build_tx_entry(best_row, 'False')  # TODO: differentiate excellent
            metrics['Alignment Modes'][aln_mode] += 1
        return best_rows is not None, gene_seen

    def find_best_score(scored_df, index_id):
        """
        Finds the best transcript in the pre-sorted filtered DataFrame, handling the case where it does not exist.
        """
        try:
            tx_rows = scored_df.xs(index_id)
            best_score = tx_rows.iloc[0].ConsensusScore
            return tx_rows[tx_rows.ConsensusScore == best_score]
        except KeyError:
            return None

    def build_tx_entry(best_row, failed_gene):
        """
        Constructs a dictionary of all of the useful attributes we have gleaned about this consensus transcript
        :param best_row:
        :return: dict of attributes
        """
        d = {'source_transcript': best_row.name, 'source_gene': best_row.Gene,
             'transcript_mode': tools.nameConversions.alignment_type(best_row.name), 'score': best_row.ConsensusScore,
             'failed_gene': failed_gene}
        return d

    def rescue_missing_gene(scored_mrna_df, scored_cds_df):
        """
        Attempts to resolve a missing gene by looking at the single best transcript within the same biotype.
        If the biotype is protein_coding, favors the CDS version
        :param scored_mrna_df: mRNA alignment mode DataFrame
        :param scored_cds_df: CDS alignment mode DataFrame
        :return: best_id or None
        """
        # combine the dataframes
        df = pd.merge(scored_mrna_df.reset_index(), scored_cds_df.reset_index(), how='outer', indicator='AlignmentMode')
        df = df.sort_values('ConsensusScore', ascending=False)
        # favor the parent biotype, where possible
        biotype_df = df[df.GeneBiotype == df.TranscriptBiotype]
        if len(biotype_df) != 0:
            df = biotype_df
        df = df.set_index('GeneId')
        best_rows = find_best_score(df, df.index[0])
        if best_rows is not None:
            metrics['Gene Rescue'] += 1
            metrics['Transcript Modes'][evaluate_ties(best_rows)] += 1
            aln_id, best_row = best_rows.iterrows().next()
            consensus[aln_id] = build_tx_entry(best_row, 'True')

    # load all genePreds
    tx_dict = tools.transcripts.load_gps(args.gp_list)
    ref_tx_dict = tools.transcripts.get_gene_pred_dict(args.annotation_gp)

    # load annotation data
    ref_df = tools.sqlInterface.load_annotation(args.ref_db_path)

    # construct a map to iterate over
    gene_transcript_map = tools.sqlInterface.get_gene_transcript_map(args.ref_db_path)

    # initialize the database session
    session = tools.sqlInterface.start_session(args.db_path)

    # did we run augustusCGP? If so, the consensus finding algorithm adds a few steps
    cgp = True if 'augCGP' in args.transcript_modes else False

    mrna_modes = list(set(args.transcript_modes) - {'augCGP'})
    cds_modes = args.transcript_modes

    # store the final consensus transcripts with their tags
    consensus = {}
    # store some metrics for plotting
    metrics = {'Alignment Modes': collections.Counter(), 'Transcript Modes': collections.Counter(),
               'Gene Failed': set(), 'Transcript Failed': set(), 'Gene Rescue': set()}

    # begin consensus finding
    for gene_id, tx_set in gene_transcript_map.iteritems():
        # extract all classification data for these transcripts
        mrna_df = load_metrics_evaluations_into_merged_table(session, mrna_modes, gene_id, ref_tx_dict, tx_dict, ref_df,
                                                             aln_mode='mRNA')
        cds_df = load_metrics_evaluations_into_merged_table(session, cds_modes, gene_id, ref_tx_dict, tx_dict, ref_df,
                                                            aln_mode='CDS')

        # score the alignments
        scored_mrna_df = score_df(mrna_df)
        scored_cds_df = score_df(cds_df)

        # filter the alignments
        filtered_mrna = filter_alignments(scored_mrna_df, aln_mode='mRNA')
        filtered_cds = filter_alignments(scored_cds_df, aln_mode='CDS')

        # indicator variable: did this gene get at least one transcript of its own biotype included?
        gene_seen = False

        # begin iterating over transcripts for this gene, trying first mRNA then CDS
        for tx_id in tx_set:
            tx_included, gene_seen = evaluate_tx_alns(filtered_mrna, tx_id, 'mRNA', gene_seen)
            if tx_included:
                tx_included, gene_seen = evaluate_tx_alns(filtered_cds, tx_id, 'CDS', gene_seen)
                if tx_included:
                    metrics['Transcript Failed'].add(tx_id)

        # attempt gene rescue, pick the one transcript to represent this
        if gene_seen is False:
            metrics['Gene Failed'].add(gene_id)
            rescue_missing_gene(mrna_df, cds_df)

        # cgp-specific novel introns
        if cgp:
            novel_isoforms, transcripts_to_include = find_novel_cgp_splices(tx_dict, consensus, cds_df,
                                                                            gene_transcript_map)
            metrics['CGP']['Novel isoforms'] = novel_isoforms
            consensus.update(transcripts_to_include)

    if cgp:
        novel_transcripts = {tx.name: tx.name for tx in tx_dict.itervalues() if 'jg' in tx.name2}
        novel_genes = {tx.name2 for tx in tx_dict.itervalues() if 'jg' in tx.name2}
        metrics['CGP'] = {'Novel genes': len(novel_genes), 'Novel transcripts': len(novel_transcripts)}
        consensus.update(novel_transcripts)

    write_consensus(args, consensus, tx_dict)
    return metrics


###
# Data loading and merging functions
###


def load_metrics_evaluations_into_merged_table(session, tx_modes, gene_id, ref_tx_dict, tx_dict, ref_df,
                                               aln_mode='mRNA'):
    """
    Loads all of the metrics and evaluations for all alignments associated with a given gene_id, bringing in outside
    data from the transcripts themselves.
    :param session:
    :param tx_modes:
    :param gene_id:
    :param ref_tx_dict:
    :param tx_dict:
    :param ref_df:
    :param aln_mode:
    :return: DataFrame
    """
    def reduce_intron_vectors(aln_ids, intron_vectors):
        """intron vector is stored as a comma separated string. Reduce this, taking aln_mode into account"""
        r = []
        for aln_id, intron_vector in zip(*[aln_ids, intron_vectors]):
            num_supported = 0
            scores = map(int, intron_vector.split(','))
            tx = tx_dict[aln_id]
            for intron, score in zip(*[tx.intron_intervals, scores]):
                if aln_mode == 'CDS' and not intron.subset(tx.coding_interval):  # don't look at this intron
                    continue
                if score == 0:  # this intron is not supported
                    continue
                num_supported += 1
            r.append(num_supported)
        return r

    def calculate_num_introns(tx):
        """calculates the number of introns we are looking at based on aln_mode"""
        if aln_mode == 'mRNA':
            return len(tx.intron_intervals)
        else:
            return len([x for x in tx.intron_intervals if x.subset(tx.coding_interval)])

    dfs = []
    intron_dfs = []
    # load the database tables for this gene
    for tx_mode in tx_modes:
        metrics_table = tools.sqlInterface.tables[aln_mode][tx_mode]['metrics']
        evaluations_table = tools.sqlInterface.tables[aln_mode][tx_mode]['evaluation']
        intron_table = tools.sqlInterface.tables['hgm'][tx_mode]
        mc_df = tools.sqlInterface.load_metrics(metrics_table, session, gene_id)
        ec_df = tools.sqlInterface.load_evaluation(evaluations_table, session, gene_id)
        intron_df = tools.sqlInterface.load_intron_vector(intron_table, session, gene_id)
        dfs.extend([mc_df, ec_df])
        intron_dfs.append(intron_df)

    # combine tables, pivot, merge
    eval_df = pd.concat(dfs)
    eval_df = pd.pivot_table(eval_df, index=['TranscriptId', 'AlignmentId'], columns='classifier', values='value',
                             fill_value=0)
    intron_df = pd.concat(intron_dfs)
    merged_df = pd.merge(eval_df.reset_index(), intron_df, on=['TranscriptId', 'AlignmentId'])
    merged_df = pd.merge(merged_df, ref_df, on=['GeneId', 'TranscriptId'], suffixes=['_Tgt', '_Ref'])
    # add in columns based on aln_mode, tx_dict, ref_tx_dict
    merged_df['NumReferenceIntrons'] = [len(ref_tx_dict[tx].intron_intervals) for tx in merged_df.TranscriptId]
    merged_df['NumIntrons'] = [calculate_num_introns(tx_dict[tx]) for tx in merged_df.AlignmentId]
    merged_df['NumSupportedIntrons'] = reduce_intron_vectors(merged_df.AlignmentId, merged_df.IntronVector)
    df = merged_df.set_index(['TranscriptId', 'AlignmentId'])
    return df


def find_novel_cgp_splices(tx_dict, consensus, cds_df, gene_transcript_map):
    """
    Finds novel splice junctions in CGP transcripts. If there are any, these get included as a novel isoform.
    :param tx_dict: Combined dictionary of GenePredTranscript objects
    :param consensus: Dict mapping {tx_id: aln_id} for transcripts that are in the consensus set
    :param cds_df:
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

    def find_existing_splices(tx_dict, consensus, tx_ids):
        existing_splices = set()
        for tx in tx_ids:
            try:  # this transcript may not be in the consensus
                consensus_tx = tx_dict[consensus[tx]]
            except KeyError:
                continue
            existing_splices.update(consensus_tx.intron_intervals)
        return existing_splices

    def find_supported_intervals(tx, cds_df):
        intervals = set()
        intron_vector = cds_df.ix[tx.name].IntronVector
        for interval, intron_score in zip(*[tx.intron_intervals, intron_vector]):
            if intron_score > 0:
                intervals.add(interval)
        return intervals

    novel_isoforms = 0
    transcripts_to_include = {}
    consensus_tx_set = set(consensus.itervalues())
    cgp_dict = create_cgp_dict(consensus_tx_set, tx_dict)
    for gene_id, tx_ids in gene_transcript_map.iteritems():
        existing_splices = find_existing_splices(tx_dict, consensus, tx_ids)
        for tx in cgp_dict[gene_id]:
            tx_supported_intervals = find_supported_intervals(tx, cds_df)
            if len(tx_supported_intervals - existing_splices) > 0:
                novel_isoforms += 1
                transcripts_to_include[tx.name] = tx.name
    return novel_isoforms, transcripts_to_include


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


def deduplicate_consensus(consensus, tx_dict, unfiltered_df):
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
    for tx_id, aln_id in consensus.iteritems():
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
    df = df.assign(ConsensusScore=scores['ConsensusScore'], StructureScore=scores['StructureScore'])
    return df.sort_values('ConsensusScore', ascending=False)


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
    structure score = 0.5 * (1 - # missing introns / # parent introns) + 0.5 * (1 - # missing exons / # parent exons)
    This function is a weighted average of the structural changes seen in the alignments.
    Parent introns/exons will be adjusted for CDS alignments to only include coding introns/exons by the sqlInterface
    query.
    Exists on the range(0, 1)

    :param s: Pandas Series
    :return: float between 0 and 1
    """
    present_introns = 1 - tools.mathOps.format_ratio(s.NumMissingIntrons, s.NumReferenceIntrons, resolve_nan=1)
    present_exons = 1 - tools.mathOps.format_ratio(s.NumMissingExons, s.NumReferenceIntrons + 1)
    score = 0.5 * present_introns + 0.5 * present_exons
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
# Outputs
###


def write_consensus(args, consensus, tx_dict):
    """
    Write the resulting gp + gp_info
    :param args:
    :param consensus:
    :return:
    """
    with args.consenus_gp.open('w') as out_gp, args.consensus_gp_info.open('w') as out_gp_info:

