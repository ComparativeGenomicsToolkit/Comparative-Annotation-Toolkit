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

Filtering transcripts:

Before scoring, the transcripts are filtered based on a series of minimum cutoffs. These are:
AlnIdentity: mRNA: >70% CDS: >80%
AlnCoverage: mRNA: >50% CDS: >90%
Intron Inequality: 2 * # missing original introns <= # reference introns - 1 or # reference introns < 5
PercentUnknownBases: mRNA: <5% CDS: <2%

These filters throw out low quality alignments via the identity/coverage filter as well as the percent unknown bases.
The intron inequality is an important filter to preventing retroposed pseudogenes from being assigned as the parent.

If the filter process throws out all transcripts for a gene, then one transcript with the best consensus score will be
used, regardless of filtering. This transcript will be marked as failing.

GFF3 tags generated in this process:
1. source_transcript: The name of the parent transcript, if it exists
2. source_gene: The name of the parent gene, if it exists
3. source_gene_common_name: The common name of the parent gene, if it is different from the source gene
4. transcript_mode: The name of the mode of operation that generated this transcript
5. score: the consensus score of the transcript
6. category: the transcript category. One of excellent, pass, fail, novel
7. alternative_source_transcripts: A comma separated list of alternate IDs for this transcript
8. failed_gene: This transcript is the single representative for a failed transcript
9. novel_sequence: This transcript is a novel gene prediction produced by CGP
10. paralogy: The number of paralogs that were mapped over when transMap mapped this transcript
11. gene_biotype: gene biotype
12. transcript_biotype: transcript biotype
"""
import copy
import collections

import pandas as pd

import tools.intervals
import tools.mathOps
import tools.fileOps
import tools.sqlInterface
import tools.transcripts
import tools.nameConversions
from tools.defaultOrderedDict import DefaultOrderedDict

cds_cutoffs = {'AlnIdentity': 0.8, 'AlnCoverage': 0.9, 'PercentUnknownBases': 0.02}
mrna_cutoffs = {'AlnIdentity': 0.7, 'AlnCoverage': 0.5, 'PercentUnknownBases': 0.05}
id_template = '{genome:.10}_{tag_type}{unique_id:07d}'


def generate_consensus(args, genome):
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
    # load all genePreds
    tx_dict = tools.transcripts.load_gps(args.gp_list)
    ref_tx_dict = tools.transcripts.get_gene_pred_dict(args.annotation_gp)

    # load annotation data
    ref_df = tools.sqlInterface.load_annotation(args.ref_db_path)
    common_name_map = dict(zip(*[ref_df.GeneId, ref_df.GeneName]))

    # load transMap results
    tm_eval = tools.sqlInterface.load_alignment_evaluation(args.db_path)
    tm_eval = tm_eval.set_index('AlignmentId')

    # construct a map to iterate over
    gene_transcript_map = tools.sqlInterface.get_gene_transcript_map(args.ref_db_path)

    # did we run augustusCGP? If so, the consensus finding algorithm adds a few steps
    cgp = True if 'augCGP' in args.transcript_modes else False

    mrna_modes = list(set(args.transcript_modes) - {'augCGP'})
    cds_modes = args.transcript_modes

    # initialize the database session
    session = tools.sqlInterface.start_session(args.db_path)
    # load all data into dataframes
    mrna_intron_df, mrna_df = load_metrics_evaluations_into_merged_table(session, mrna_modes, ref_tx_dict,
                                                                         tx_dict, ref_df, aln_mode='mRNA')
    cds_intron_df, cds_df = load_metrics_evaluations_into_merged_table(session, cds_modes, ref_tx_dict,
                                                                       tx_dict, ref_df, aln_mode='CDS')

    # store the final consensus transcripts with their tags
    consensus = {}
    # store some metrics for plotting
    metrics = {'Alignment Modes': collections.Counter(), 'Transcript Modes': collections.Counter(),
               'Gene Failed': 0, 'Transcript Failed': 0, 'Gene Rescue': 0, 'Gene Missing': 0,
               'Duplicate transcripts': 0}

    # novel genes/transcripts from CGP based on parent assignment
    if cgp:
        find_novel_transcripts(tx_dict, consensus, metrics)

    # begin consensus finding
    for gene_id, tx_set in gene_transcript_map.iteritems():
        # slice out the classification data for this gene, handling genes that never mapped
        try:
            gene_mrna_df = mrna_df.xs(gene_id)
        except KeyError:
            gene_mrna_df = pd.DataFrame()
        try:
            gene_cds_df = cds_df.xs(gene_id)
        except KeyError:
            gene_cds_df = pd.DataFrame()

        if len(gene_cds_df) == len(gene_mrna_df) == 0:
            metrics['Gene Missing'] += 1
            continue

        # score the alignments
        scored_mrna_df = score_df(gene_mrna_df)
        scored_cds_df = score_df(gene_cds_df)

        # filter the alignments
        filtered_mrna = filter_alignments(scored_mrna_df, aln_mode='mRNA')
        filtered_cds = filter_alignments(scored_cds_df, aln_mode='CDS')

        # indicator variable: did this gene get at least one transcript of its own biotype included?
        gene_seen = False
        # if a tx gets included, store its ID here so that we can evaluate CGP for novel isoforms
        gene_tx_objs = []

        # begin iterating over transcripts for this gene, trying first mRNA then CDS
        for tx_id in tx_set:
            best_rows = find_best_score(filtered_mrna, tx_id)
            if best_rows is not None:
                aln_id, gene_seen = incorporate_tx(best_rows, consensus, metrics, gene_seen, 'mRNA', tx_id, gene_id,
                                                   tm_eval)
                gene_tx_objs.append(tx_dict[aln_id])
            else:  # try CDS
                best_rows = find_best_score(filtered_cds, tx_id)
                if best_rows is not None:
                    aln_id, gene_seen = incorporate_tx(best_rows, consensus, metrics, gene_seen, 'CDS', tx_id, gene_id,
                                                       tm_eval)
                    gene_tx_objs.append(tx_dict[aln_id])
                else:
                    metrics['Transcript Failed'] += 1

        # attempt gene rescue, pick the one transcript to represent this
        if gene_seen is False and (len(scored_mrna_df) > 0 or len(scored_cds_df) > 0):
            metrics['Gene Failed'] += 1
            rescue_missing_gene(scored_mrna_df, scored_cds_df, metrics, tm_eval, consensus, gene_id)

        # cgp-specific novel introns
        if cgp:
            find_novel_cgp_splices(gene_tx_objs, tx_dict, cds_intron_df.xs(gene_id), gene_id, common_name_map,
                                   consensus, gene_seen, metrics)

    # perform final filtering steps
    deduplicated_consensus = deduplicate_consensus(consensus, tx_dict, metrics)
    deduplicated_strand_resolved_consensus = resolve_opposite_strand(deduplicated_consensus, tx_dict)

    # write out reuslts. consensus tx dict has the unique names
    consensus_gene_dict = write_consensus_gps(args.consensus_gp, args.consensus_gp_info,
                                              deduplicated_strand_resolved_consensus, tx_dict, genome)
    write_consensus_gff3(consensus_gene_dict, args.consensus_gff3)

    return metrics


###
# Consensus finding functions
###


def find_novel_transcripts(tx_dict, consensus, metrics):
    """Finds novel transcripts, builds their attributes"""
    novel_transcripts = {tx.name: tx.name for tx in tx_dict.itervalues() if 'jg' in tx.name2}
    novel_genes = {tx.name2 for tx in tx_dict.itervalues() if 'jg' in tx.name2}
    metrics['CGP'] = {'Novel genes': len(novel_genes), 'Novel transcripts': len(novel_transcripts)}
    for novel_tx in novel_transcripts:
        consensus[novel_tx] = {'category': 'novel', 'novel_sequence': True, 'gene_biotype': 'unknown_likely_coding',
                               'transcript_biotype': 'unknown_likely_coding'}


def find_best_score(scored_df, index_id):
    """
    Finds the best transcript in the pre-sorted filtered DataFrame, handling the case where it does not exist.
    """
    try:
        tx_rows = scored_df.xs(index_id)
    except KeyError:
        return None
    try:
        if isinstance(tx_rows, pd.core.series.Series):
            return pd.DataFrame([tx_rows])
        else:
            best_score = tx_rows.iloc[0].ConsensusScore
    except IndexError:
        return None
    return tx_rows[tx_rows.ConsensusScore == best_score]


def incorporate_tx(best_rows, consensus, metrics, gene_seen, aln_mode, tx_id, gene_id, tm_eval):
    """incorporate a transcript into the consensus set, storing metrics. Updates gene_seen."""
    metrics['Transcript Modes'][evaluate_ties(best_rows)] += 1
    _, best_series = best_rows.iterrows().next()
    aln_id = best_series.AlignmentId
    consensus[aln_id] = build_tx_entry(best_series, False, tx_id, gene_id, tm_eval)
    metrics['Alignment Modes'][aln_mode] += 1
    return aln_id, evaluate_seen(gene_seen, best_series)


def evaluate_seen(gene_seen, best_row):
    """a gene is only seen if the transcript biotype matches the gene biotype"""
    if gene_seen is False:
        gene_seen = True if best_row.TranscriptBiotype == best_row.GeneBiotype else False
    return gene_seen


def evaluate_ties(best_rows):
    """Find out how many transcript modes agreed on this"""
    return ','.join(sorted(set([tools.nameConversions.alignment_type(x) for x in best_rows.AlignmentId])))


def build_tx_entry(best_series, failed_gene, tx_id, gene_id, tm_eval):
    """
    Constructs a dictionary of all of the useful attributes we have gleaned about this consensus transcript
    """
    d = {'source_transcript': tx_id,
         'source_gene': gene_id,
         'transcript_mode': tools.nameConversions.alignment_type(best_series.AlignmentId),
         'score': round(best_series.ConsensusScore, 2),
         'failed_gene': failed_gene,
         'category': best_series.Category,
         'gene_biotype': best_series.GeneBiotype,
         'transcript_biotype': best_series.TranscriptBiotype}
    paralogy = tm_eval.ix[tools.nameConversions.remove_augustus_alignment_number(best_series.AlignmentId)].Paralogy
    if paralogy > 0:
        d['paralogy'] = paralogy
    if best_series.GeneName != gene_id:
        d['source_gene_common_name'] = best_series.GeneName
    return d


def rescue_missing_gene(scored_mrna_df, scored_cds_df, metrics, tm_eval, consensus, gene_id):
    """
    Attempts to resolve a missing gene by looking at the single best transcript within the same biotype.
    """
    # combine the dataframes
    df = pd.merge(scored_mrna_df.reset_index(), scored_cds_df.reset_index(), how='outer', indicator='AlignmentMode',
                  copy=False)
    df = df.sort_values('ConsensusScore', ascending=False)
    # only pick if the parent biotype matches
    biotype_df = df[df.GeneBiotype == df.TranscriptBiotype]
    if len(biotype_df) == 0:
        return
    # below code is a hack to be able to make use of find_best_score()
    biotype_df['GeneId'] = [gene_id] * len(biotype_df)
    biotype_df = biotype_df.set_index('GeneId')
    best_rows = find_best_score(biotype_df, biotype_df.index[0])
    if best_rows is not None:
        metrics['Gene Rescue'] += 1
        metrics['Transcript Modes'][evaluate_ties(best_rows)] += 1
        _, best_series = best_rows.iterrows().next()
        consensus[best_series.AlignmentId] = build_tx_entry(best_series, True, best_series.TranscriptId, gene_id,
                                                            tm_eval)


def find_novel_cgp_splices(gene_tx_objs, tx_dict, cds_intron_df, gene_id, common_name_map, consensus, gene_seen,
                           metrics):
    """
    Finds novel splice junctions in CGP transcripts. If there are any, these get included as a novel isoform.
    """
    metrics['Novel isoforms'] = 0
    existing_splices = set()
    for consensus_tx in gene_tx_objs:
        existing_splices.update(consensus_tx.intron_intervals)
    cgp_txs = {tx_id for tx_id in cds_intron_df.index.get_level_values('TranscriptId')
               if tools.nameConversions.aln_id_is_cgp(tx_id)}
    intron_df = cds_intron_df.ix[cgp_txs]
    for cgp_tx in cgp_txs:
        cgp_tx_obj = tx_dict[cgp_tx]
        intron_series = intron_df.ix[cgp_tx]
        intron_vector = map(int, intron_series.IntronVector.split(','))
        for interval, intron_score in zip(*[cgp_tx_obj.intron_intervals, intron_vector]):
            if intron_score > 0 and interval not in existing_splices:
                metrics['Novel isoforms'] += 1
                consensus[cgp_tx] = {'category': 'novel', 'source_gene': gene_id,
                                     'failed_gene': not gene_seen, 'transcript_mode': 'augCGP',
                                     'transcript_biotype': 'unknown_likely_coding',
                                     'gene_biotype': 'unknown_likely_coding'}
                common_name = common_name_map[gene_id]
                if common_name != gene_id:
                    consensus[cgp_tx]['source_gene_common_name'] = common_name


###
# Data loading and merging functions
###


def load_metrics_evaluations_into_merged_table(session, tx_modes, ref_tx_dict, tx_dict, ref_df, aln_mode):
    """
    Loads all of the metrics and evaluations for all alignments associated with a given gene_id, bringing in outside
    data from the transcripts themselves.
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
        mc_df = tools.sqlInterface.load_metrics(metrics_table, session)
        ec_df = tools.sqlInterface.load_evaluation(evaluations_table, session)
        intron_df = tools.sqlInterface.load_intron_vector(intron_table, session)
        dfs.extend([mc_df, ec_df])
        intron_dfs.append(intron_df)

    # combine tables
    eval_df = pd.concat(dfs)
    intron_df = pd.concat(intron_dfs)

    # pivot tables, merge
    pivot_df = pd.pivot_table(eval_df, index=['GeneId', 'TranscriptId', 'AlignmentId'], columns='classifier',
                              values='value', fill_value=0)
    merged_df = pd.merge(pivot_df.reset_index(), intron_df, on=['GeneId', 'TranscriptId', 'AlignmentId'])

    # add in columns based on aln_mode, tx_dict, ref_tx_dict
    merged_df['NumReferenceIntrons'] = [len(ref_tx_dict[tx].intron_intervals) for tx in merged_df.TranscriptId]
    merged_df['NumIntrons'] = [calculate_num_introns(tx_dict[tx]) for tx in merged_df.AlignmentId]
    merged_df['NumSupportedIntrons'] = reduce_intron_vectors(merged_df.AlignmentId, merged_df.IntronVector)
    df = pd.merge(merged_df, ref_df, on=['GeneId', 'TranscriptId'], suffixes=['_Tgt', '_Ref'], how='left')
    df = df.set_index(['GeneId', 'TranscriptId'])
    intron_df = intron_df.set_index(['GeneId', 'TranscriptId'])
    return intron_df, df


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
    if len(merged_df) == 0:
        return merged_df
    cutoffs = cds_cutoffs if aln_mode == 'CDS' else mrna_cutoffs
    return merged_df[(merged_df['AlnIdentity'] > cutoffs['AlnIdentity']) &
                     (merged_df['AlnCoverage'] > cutoffs['AlnCoverage']) &
                     (merged_df['PercentUnknownBases'] <= cutoffs['PercentUnknownBases']) &
                     ((2 * merged_df['NumMissingIntrons'] <= merged_df['NumReferenceIntrons'] - 1) |
                      (merged_df['NumReferenceIntrons'] < 5))]


def deduplicate_consensus(consensus, tx_dict, metrics):
    """
    In the process of consensus building, we may find that we have ended up with more than one transcript for a gene
    that are actually identical. Remove these, picking the best based on their score, favoring the transcript
    whose biotype matches the parent.
    """
    def resolve_duplicate(tx_list, consensus):
        biotype_txs = [tx for tx in tx_list if
                       consensus[tx].get('gene_biotype', None) == consensus[tx].get('transcript_biotype', None)]
        if len(biotype_txs) > 0:
            sorted_scores = sorted([[tx, consensus[tx].get('score', 0)]for tx in biotype_txs], key=lambda (tx, s): -s)
            return sorted_scores[0][0]
        else:
            sorted_scores = sorted([[tx, consensus[tx].get('score', 0)]for tx in tx_list], key=lambda (tx, s): -s)
            return sorted_scores[0][0]

    def add_duplicate_field(best_tx, tx_list, consensus, deduplicated_consensus):
        deduplicated_consensus[best_tx] = consensus[best_tx]
        deduplicated_consensus[best_tx]['alternative_source_transcripts'] = ','.join(set(tx_list) - {best_tx})

    # build a dictionary mapping duplicates making use of hashing intervals
    duplicates = collections.defaultdict(list)
    for aln_id in consensus:
        tx = tx_dict[aln_id]
        duplicates[frozenset(tx.exon_intervals)].append(aln_id)

    # begin iterating
    deduplicated_consensus = {}
    for tx_list in duplicates.itervalues():
        if len(tx_list) > 1:
            metrics['Duplicate transcripts'] += 1
            best_tx = resolve_duplicate(tx_list, consensus)
            add_duplicate_field(best_tx, tx_list, consensus, deduplicated_consensus)
        else:
            tx_id = tx_list[0]
            deduplicated_consensus[tx_id] = consensus[tx_id]

    # sort by genomic interval for prettily increasing numbers
    deduplicated_consensus = sorted(deduplicated_consensus.iteritems(),
                                    key=lambda (tx, attrs): (tx_dict[tx].chromosome, tx_dict[tx].start))
    return deduplicated_consensus


def resolve_opposite_strand(deduplicated_consensus, tx_dict):
    """
    Resolves situations where multiple transcripts of the same gene are on opposite strands. Does so by looking for
    the largest sum of scores.
    """
    gene_dict = collections.defaultdict(list)
    for tx_id, attrs in deduplicated_consensus:
        tx_obj = tx_dict[tx_id]
        gene_dict[tx_obj.name2].append([tx_obj, attrs])

    deduplicated_strand_resolved_consensus = []
    for gene in gene_dict:
        tx_objs, attrs = zip(*gene_dict[gene])
        if len(set(tx_obj.strand for tx_obj in tx_objs)) > 1:
            strand_scores = collections.Counter()
            for tx_obj, attrs in gene_dict[gene]:
                strand_scores[tx_obj.strand] += attrs['score']
            best_strand = sorted(strand_scores.items())[0][0]
            for tx_obj, attrs in gene_dict[gene]:
                if tx_obj.strand == best_strand:
                    deduplicated_strand_resolved_consensus.append([tx_obj.name, attrs])
        else:
            deduplicated_strand_resolved_consensus.extend([[tx_obj.name, attrs] for tx_obj, attrs in gene_dict[gene]])
    return deduplicated_strand_resolved_consensus


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
        consensus_score, percent_junction_support = calculate_consensus_score(s)
        scores['ConsensusScore'].append(consensus_score)
        scores['Category'].append(categorize_aln(s, percent_junction_support))
    df = df.assign(ConsensusScore=scores['ConsensusScore'], Category=scores['Category'])
    return df.sort_values('ConsensusScore', ascending=False)


def categorize_aln(s, percent_junction_support):
    """
    Categorizes an alignment as excellent/pass/fail based on their structure.
    A excellent alignment has:
    1) No coding indels, 2) No in frame stops, 3) Complete ends (unless incomplete in the reference), 4) all junctions
    supported.
    A passing alignment has:
    1) No in frame stops, 2) 75% of junctions supported
    :param s: pandas Series
    :param percent_junction_support: percentage as a float
    :return: string
    """
    if s.TranscriptBiotype == 'protein_coding':
        if (s.StartCodon_Tgt == 1 or s.StartCodon_Ref == 0) and (s.StopCodon_Tgt == 1 or s.StopCodon_Ref == 0) and \
              s.get('InFrameStop', 0) == 0 and s.get('CodingDeletion', 0) == 0 and s.get('CodingInsertion', 0) == 0 \
              and percent_junction_support == 1:
            return 'Excellent'
        elif s.get('InFrameStop', 0) == 0 and percent_junction_support > 0.75:
            return 'Passing'
    else:
        if percent_junction_support == 1:
            return 'Excellent'
        if percent_junction_support > 0.75:
            return 'Passing'
    return 'Fail'


def calculate_consensus_score(s):
    """
    consensus score = 50 * (1 - badness) + 25 * structure score + 25 * # supported junctions
    Exists in the range (0, 100). Higher is better.

    coding consensus score = consensus score + cds bonus + evaluation bonus - evaluation penalties
    For coding transcripts consensus.py modifies the consensus score by removing as much as 5 and adding as much as 4

    Exists in the range (-5, 104)

    :param s: Pandas Series
    :return: float
    """
    percent_junction_support = supported_junctions(s)
    score = 50 * (1 - s.Badness) + 25 * calculate_structure_score(s) + 25 * percent_junction_support
    if s.TranscriptBiotype == 'protein_coding':
        score = score + evaluation_bonus(s) - evaluation_penalties(s)
        assert -5 <= score <= 104
    else:
        assert 0 <= score <= 100
    return score, percent_junction_support


def indicator(item):
    """Pass a Series member and returns a 1 if the value in the member is > 0"""
    return 1 if item > 0 else 0


def calculate_structure_score(s):
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


def write_consensus_gps(consensus_gp, consensus_gp_info, deduplicated_strand_resolved_consensus, tx_dict, genome):
    """
    Write the resulting gp + gp_info, generating genome-specific unique identifiers
    """
    gene_count = 0
    tx_count = 1
    consensus_gene_dict = DefaultOrderedDict(lambda: DefaultOrderedDict(list))  # used to make gff3 next
    gp_infos = []
    genes_seen = set()
    with consensus_gp.open('w') as out_gp:
        for tx, attrs in deduplicated_strand_resolved_consensus:
            tx_obj = copy.deepcopy(tx_dict[tx])
            tx_obj.name = id_template.format(genome=genome, tag_type='T', unique_id=tx_count)
            tx_count += 1
            if tx_obj.name2 not in genes_seen:
                genes_seen.add(tx_obj.name2)
                gene_count += 1
            tx_obj.name2 = id_template.format(genome=genome, tag_type='G', unique_id=gene_count)
            out_gp.write('\t'.join(tx_obj.get_gene_pred()) + '\n')
            consensus_gene_dict[tx_obj.chromosome][tx_obj.name2].append([tx_obj, attrs])
            gp_info = attrs.copy()
            gp_info['transcript_id'] = tx_obj.name
            gp_info['gene_id'] = tx_obj.name2
            gp_infos.append(gp_info)
    gp_info_df = pd.DataFrame(gp_infos)
    gp_info_df = gp_info_df.set_index(['gene_id', 'transcript_id'])
    gp_info_df.to_csv(consensus_gp_info.path, sep='\t')
    return consensus_gene_dict


def write_consensus_gff3(consensus_gene_dict, consensus_gff3):
    """
    Write the consensus set in gff3 format
    """
    def convert_frame(exon_frame):
        """converts genePred-style exonFrame to GFF-style phase"""
        mapping = {0: 0, 1: 2, 2: 1}
        return mapping[exon_frame]

    def convert_attrs(attrs, id_field):
        """converts the attrs dict to a attributes field. assigns name to the gene common name for display"""
        attrs['ID'] = id_field
        if 'source_gene_common_name' in attrs:
            attrs['Name'] = attrs['source_gene_common_name']
        attrs_str = ['='.join([key, str(val)]) for key, val in sorted(attrs.iteritems())]
        return ';'.join(attrs_str)

    def generate_gene_record(chrom, tx_objs, gene_id, attrs):
        """calculates the gene interval for this list of tx"""
        intervals = set()
        for tx in tx_objs:
            intervals.update(tx.exon_intervals)
        intervals = sorted(intervals)
        strand = tx_objs[0].strand
        # subset the attrs to gene fields
        useful_keys = ['source_gene_common_name', 'source_gene', 'gene_biotype', 'failed_gene']
        attrs = {key: attrs[key] for key in useful_keys if key in attrs}
        attrs_field = convert_attrs(attrs, gene_id)
        return [chrom, 'CAT', 'gene', intervals[0].start + 1, intervals[-1].stop + 1, '.', strand, '.', attrs_field]

    def generate_transcript_record(chrom, tx_obj, attrs):
        """generates transcript records, calls generate_exon_records to generate those too"""
        tx_id = tx_obj.name
        gene_id = tx_obj.name2
        attrs['Parent'] = gene_id
        attrs_field = convert_attrs(attrs, tx_id)
        yield [chrom, 'CAT', 'transcript', tx_obj.start + 1, tx_obj.stop + 1, '.', tx_obj.strand, '.', attrs_field]
        for line in generate_exon_records(chrom, tx_obj, tx_id, attrs):
            yield line
        for line in generate_start_stop_codon_records(chrom, tx_obj, tx_id, attrs):
            yield line

    def generate_exon_records(chrom, tx_obj, tx_id, attrs):
        """generates exon records"""
        attrs['Parent'] = tx_id
        for i, (exon, exon_frame) in enumerate(zip(*[tx_obj.exon_intervals, tx_obj.exon_frames]), 1):
            attrs_field = convert_attrs(attrs, 'exon:{}:{}'.format(tx_id, i))
            yield [chrom, 'CAT', 'exon', exon.start + 1, exon.stop + 1, '.', exon.strand, '.', attrs_field]
            cds_interval = exon.intersection(tx_obj.coding_interval)
            if cds_interval is not None:
                attrs_field = convert_attrs(attrs, 'CDS:{}:{}'.format(tx_id, i))
                yield [chrom, 'CAT', 'CDS', cds_interval.start + 1, cds_interval.stop + 1, '.', exon.strand,
                       convert_frame(exon_frame), attrs_field]

    def generate_start_stop_codon_records(chrom, tx_obj, tx_id, attrs):
        """generate start/stop codon GFF3 records, handling frame appropriately"""
        cds_frames = [x for x in tx_obj.exon_frames if x != -1]
        if tx_obj.cds_start_stat == 'cmpl':
            attrs_field = convert_attrs(attrs, 'start_codon:{}'.format(tx_id))
            start, stop = tools.transcripts.get_start_interval(tx_obj)
            if tx_obj.strand == '-':
                start_frame = convert_frame(cds_frames[-1])
            else:
                start_frame = convert_frame(cds_frames[0])
            yield [chrom, 'CAT', 'start_codon', start + 1, stop + 1, '.', tx_obj.strand, start_frame, attrs_field]
        if tx_obj.cds_end_stat == 'cmpl':
            attrs_field = convert_attrs(attrs, 'stop_codon:{}'.format(tx_id))
            start, stop = tools.transcripts.get_stop_interval(tx_obj)
            if tx_obj.strand == '-':
                stop_frame = convert_frame(cds_frames[-1])
            else:
                stop_frame = convert_frame(cds_frames[0])
            yield [chrom, 'CAT', 'stop_codon', start + 1, stop + 1, '.', tx_obj.strand, stop_frame, attrs_field]

    # main gff3 writing logic
    with consensus_gff3.open('w') as out_gff3:
        out_gff3.write('##gff-version 3\n')
        for chrom in consensus_gene_dict:
            for gene_id, tx_list in consensus_gene_dict[chrom].iteritems():
                tx_objs, attrs_list = zip(*tx_list)
                attrs = tx_list[0][1]  # grab the attrs from the first transcript
                tools.fileOps.print_row(out_gff3, generate_gene_record(chrom, tx_objs, gene_id, attrs))
                tx_lines = []
                for tx_obj, attrs in tx_list:
                    tx_lines.extend(list(generate_transcript_record(chrom, tx_obj, attrs)))
                tx_lines = sorted(tx_lines, key=lambda l: l[3])
                tools.fileOps.print_rows(out_gff3, tx_lines)
