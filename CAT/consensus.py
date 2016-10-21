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
6. alternative_source_transcripts: A comma separated list of alternate IDs for this transcript
7. failed_gene: This transcript is the single representative for a failed transcript
8. transcript_class: One of failed, passing, excellent, novel
9. paralogy: The number of paralogs that were mapped over when transMap mapped this transcript
10. gene_biotype: gene biotype
11. transcript_biotype: transcript biotype
12. paralog_status: confident if the paralog was confidently resolved, not_confident if it was not
13. alternative_source_transcripts: Other possible transcripts, if this was collapsed as the result of deduplication
14. gene_alternate_contigs: If the --resolve-split-genes flag was set, contigs that this gene was also found on are
    comma separated in this tag.
"""
import collections
import copy
import luigi
import math

import pandas as pd
import numpy as np

import tools.intervals
import tools.mathOps
import tools.fileOps
import tools.sqlInterface
import tools.transcripts
import tools.nameConversions
from tools.defaultOrderedDict import DefaultOrderedDict

pd.options.mode.chained_assignment = None  # false positive in novel splice finding code
id_template = '{genome:.10}_{tag_type}{unique_id:07d}'


def generate_consensus(args, genome):
    """
    Entry point for consensus finding algorithm. Main consensus finding logic is here.

    A mega DataFrame is constructed with the metrics for all transcript alignments in both coding and non-coding.
    Coding alignments are transcript-transcript BLAT for all transcript modes and alignment modes (mRNA and CDS).
    Noncoding alignments are simply transMap output. See their individual scoring functions for more details.

    The first step is to incorporate any CGP transcripts that were not assigned a parental gene as novel genes.
    These are filtered for a minimum number of exons and intron support.

    The scored dataframe is sorted first alignment mode then by consensus score. This means that we only look at the CDS
    alignments if none of the mRNA alignments are of sufficient quality. The consensus logic is as follows:
    For each gene, see if we have it. If not, it is missing. If we do have it, see if all transcripts are failing.
    If all transcripts are failing, pick one. Otherwise, take all transcripts for this gene that are not failing.

    :param args: Argument namespace from luigi
    """
    # stores a mapping of alignment IDs to tags for the final consensus set
    consensus_dict = {}

    # store some metrics for plotting
    metrics = {'Alignment Modes': collections.Counter(), 'Transcript Modes': collections.Counter(),  # coding only
               'Gene Failed': collections.Counter(), 'Transcript Failed': collections.Counter(),
               'Transcript Missing': collections.Counter(),
               'Gene Rescue': collections.Counter(), 'Gene Missing': collections.Counter(),
               'Duplicate transcripts': collections.Counter(),
               'Discarded by strand resolution': 0,
               'Novel isoforms': 0, 'Novel genes': 0,
               'Transcript Categories': collections.defaultdict(lambda: collections.Counter()),
               'Coverage': collections.defaultdict(list),
               'Identity': collections.defaultdict(list),
               'Consensus Score': collections.defaultdict(list),
               'Splice Support': collections.defaultdict(list),
               'transMap Splice Support': collections.defaultdict(list)}

    # load all genePreds
    tx_dict = tools.transcripts.load_gps(args.gp_list)
    # load annotation data
    ref_df = tools.sqlInterface.load_annotation(args.ref_db_path)

    # load transMap evaluation data
    tm_eval = load_transmap_evals(args.db_path, ref_df)
    coding_cutoff = tools.sqlInterface.load_tm_fit(args.db_path)

    # did we run augustusTMR? In other words, do we have RNA-seq hints? If so, the scoring metric changes to
    # account for the added information from Hgm
    if args.hints_db_has_rnaseq is True:
        intron_df = load_intron_vectors(args.db_path, args.transcript_modes, tx_dict, ref_df)
    else:
        intron_df = None

    scored_df = score_alignments(args.db_path, args.transcript_modes, ref_df, tm_eval, coding_cutoff, intron_df)

    # did we run augustusCGP? If so, find novel genes within cutoffs
    # TODO: expose these cutoffs to the user as options
    if args.augustus_cgp is True:
        consensus_dict.update(find_novel_transcripts(intron_df, tx_dict, metrics))
        metrics['Novel genes'] = len(consensus_dict)

    # gene transcript map to iterate over so that we capture missing gene information
    gene_transcript_map = tools.sqlInterface.get_gene_transcript_map(args.ref_db_path)
    gene_biotype_map = tools.sqlInterface.get_gene_biotype_map(args.ref_db_path)
    transcript_biotype_map = tools.sqlInterface.get_transcript_biotype_map(args.ref_db_path)
    common_name_map = dict(zip(*[ref_df.GeneId, ref_df.GeneName]))

    for gene_id, tx_list in gene_transcript_map.iteritems():
        gene_consensus_dict = {}
        gene_df = slice_df(scored_df, gene_id)
        gene_biotype = gene_biotype_map[gene_id]
        if len(gene_df) == 0:
            metrics['Gene Missing'][gene_biotype] += 1
            continue
        failed_gene = is_failed_df(gene_df)
        if failed_gene is True:
            aln_id, d = rescue_failed_gene(gene_df, tx_dict, gene_id, metrics)
            gene_consensus_dict[aln_id] = d
            metrics['Gene Failed'][gene_biotype] += 1
        else:
            for tx_id in tx_list:
                tx_biotype = transcript_biotype_map[tx_id]
                tx_df = slice_df(gene_df, tx_id)
                if len(tx_df) == 0:
                    metrics['Transcript Missing'][tx_biotype] += 1
                elif is_failed_df(tx_df):
                    metrics['Transcript Failed'][tx_biotype] += 1
                else:
                    best_rows = find_best_score(tx_df)
                    aln_id, d = incorporate_tx(best_rows, gene_id, metrics, failed_gene=True)
                    gene_consensus_dict[aln_id] = d
        if args.augustus_cgp is True:
            gene_consensus_dict.update(find_novel_cgp_splices(gene_consensus_dict, gene_df, tx_dict, gene_id,
                                                              common_name_map, metrics, failed_gene))
        consensus_dict.update(gene_consensus_dict)

    # perform final filtering steps
    deduplicated_consensus = deduplicate_consensus(consensus_dict, tx_dict, metrics)
    deduplicated_strand_resolved_consensus = resolve_opposite_strand(deduplicated_consensus, tx_dict, metrics)

    # sort by genomic interval for prettily increasing numbers
    final_consensus = sorted(deduplicated_strand_resolved_consensus,
                             key=lambda (tx, attrs): (tx_dict[tx].chromosome, tx_dict[tx].start))

    # calculate final gene set completeness
    calculate_completeness(final_consensus, metrics)

    # write out results. consensus tx dict has the unique names
    consensus_gene_dict = write_consensus_gps(args.consensus_gp, args.consensus_gp_info,
                                              final_consensus, tx_dict, genome)
    write_consensus_gff3(consensus_gene_dict, args.consensus_gff3)

    return metrics


def load_transmap_evals(db_path, ref_df):
    """
    Loads the database tables associated with both transMap filtering and transMap evaluation, merging them
    """
    # load transMap results
    tm_eval = tools.sqlInterface.load_alignment_evaluation(db_path)
    # load transMap filtering results
    tm_filter_eval = tools.sqlInterface.load_filter_evaluation(db_path)

    # combine transMap evaluation and transMap filtering into one table
    # the transMap filtering columns are used for tags in the output
    merged = pd.merge(tm_eval, tm_filter_eval, on=['TranscriptId', 'AlignmentId'])
    ref_merged = pd.merge(merged, ref_df, on=['GeneId', 'TranscriptId'])
    return ref_merged


def load_intron_vectors(db_path, tx_modes, tx_dict, ref_df):
    """
    Loads the intron vector table output by the homGeneMapping module. Returns a DataFrame with
    """
    def split_intron_string(s):
        """Intron vector is stored as a comma separated string. Turn this into a list."""
        return map(int, s.IntronVector.split(','))

    def reduce_intron_vectors(s, coding):
        """Reduce the intron vector list, dealing with coding transcripts"""
        if coding is False:
            num_supported = len([x for x in s.IntronVector if x > 0])
            return tools.mathOps.format_ratio(num_supported, s.NumIntrons)
        else:
            num_supported = 0
            tx = tx_dict[s.AlignmentId]
            for intron, score in zip(*[tx.intron_intervals, s.IntronVector]):
                if intron.subset(tx.coding_interval) and score > 0:
                    num_supported += 1
            # we give transcripts with no introns a score of 1 here
            return tools.mathOps.format_ratio(num_supported, s.NumCodingIntrons, resolve_nan=1)

    def calculate_coding_introns(s):
        """calculates the number of coding introns, if this transcript is coding"""
        if s.TranscriptBiotype != 'protein_coding':
            return None
        else:
            tx = tx_dict[s.AlignmentId]
            return len([x for x in tx.intron_intervals if x.subset(tx.coding_interval)])

    session = tools.sqlInterface.start_session(db_path)
    intron_dfs = []
    # load the database tables for this gene
    for tx_mode in tx_modes:
        intron_table = tools.sqlInterface.tables['hgm'][tx_mode]
        intron_df = tools.sqlInterface.load_intron_vector(intron_table, session)
        intron_dfs.append(intron_df)

    intron_df = pd.concat(intron_dfs)
    intron_df = pd.merge(intron_df, ref_df, on=['GeneId', 'TranscriptId'], how='left')
    # start calculating support levels for consensus finding
    intron_df['IntronVector'] = intron_df.apply(split_intron_string, axis=1)
    intron_df['NumIntrons'] = [len(tx_dict[tx].intron_intervals) for tx in intron_df.AlignmentId]
    intron_df['NumCodingIntrons'] = intron_df.apply(calculate_coding_introns, axis=1)
    intron_df['PercentIntronsSupported'] = intron_df.apply(reduce_intron_vectors, coding=False, axis=1)
    intron_df['PercentCodingIntronsSupported'] = intron_df.apply(reduce_intron_vectors, coding=True, axis=1)
    # we don't carry along transcript ID because it will conflict with CGP transcript IDs
    return intron_df[['AlignmentId', 'GeneId', 'IntronVector',
                      'PercentIntronsSupported', 'PercentCodingIntronsSupported', 'NumIntrons']]


def score_alignments(db_path, transcript_modes, ref_df, tm_eval, coding_cutoff, intron_df):
    """
    Returns a combined sorted dataframe of scored alignments.
    """
    coding_df = load_metrics_evaluations(db_path, transcript_modes, ref_df, tm_eval, coding_cutoff, intron_df)
    non_coding_df = score_non_coding_df(tm_eval, intron_df)
    combined_df = pd.concat([coding_df, non_coding_df])

    # convert the AlnMode column to be categorical for later sorting
    combined_df.AlnMode = pd.Categorical(combined_df.AlnMode, ['mRNA', 'CDS'])
    sorted_df = combined_df.sort_values(['GeneId', 'TranscriptId', 'AlnMode', 'ConsensusScore'],
                                        ascending=[True, True, True, False])  # put consensus in reverse order
    return sorted_df.set_index(['GeneId', 'TranscriptId'])


def load_metrics_evaluations(db_path, transcript_modes, ref_df, tm_eval, coding_cutoff, intron_df=None):
    """
    Loads all of the metrics and evaluations produced by transcript alignments for transMap (coding only),
    AugustusTM(R) and AugustusCGP (CDS alignment mode only).
    """
    session = tools.sqlInterface.start_session(db_path)
    dfs = []
    # load the database tables for this gene
    for aln_mode in ['mRNA', 'CDS']:
        for tx_mode in transcript_modes:
            if tx_mode == 'augCGP' and aln_mode == 'mRNA':
                continue
            metrics_table = tools.sqlInterface.tables[aln_mode][tx_mode]['metrics']
            evaluations_table = tools.sqlInterface.tables[aln_mode][tx_mode]['evaluation']
            mc_df = tools.sqlInterface.load_metrics(metrics_table, session)
            mc_df = pd.pivot_table(mc_df, index=['GeneId', 'TranscriptId', 'AlignmentId'], columns='classifier',
                                   values='value', fill_value=None).reset_index()
            mc_df['AlnMode'] = [aln_mode] * len(mc_df)
            ec_df = tools.sqlInterface.load_evaluation(evaluations_table, session)
            ec_df = pd.pivot_table(ec_df, index=['GeneId', 'TranscriptId', 'AlignmentId'], columns='classifier',
                                   values='value', fill_value=None).reset_index()
            ec_df['AlnMode'] = [aln_mode] * len(ec_df)
            dfs.extend([mc_df, ec_df])

    session.close()
    df = pd.concat(dfs)

    # coverage filter, much higher than transMap because this is for consensus
    filtered_df = df[(df.AlnCoverage > 0.4) & (df.AlnIdentity >= coding_cutoff)]

    # bring in biotype and gene information from the ref database
    ref_merged = pd.merge(filtered_df, ref_df, on=['GeneId', 'TranscriptId'], how='left', suffixes=['_Tgt', '_Ref'])

    # slice out the tm_eval columns we need to make this easier
    tm_eval_subset = tm_eval[['TranscriptId', 'ParalogStatus', 'GeneAlternateContigs', 'TranscriptClass', 'Paralogy']]
    # we use a left outer join because of CGP transcripts
    merged_df = pd.merge(ref_merged, tm_eval_subset, on='TranscriptId', how='left')

    # bring in intron support, if we have it
    if intron_df is not None:
        # we use a left outer join because of single exon transcripts
        merged_df = pd.merge(merged_df, intron_df, on=['GeneId', 'AlignmentId'], how='left')

    # score the dataframe
    scored_df = score_coding_df(merged_df, intron_df is not None)
    return scored_df


def score_coding_df(merged_df, has_rnaseq_data):
    """
    Scores a DataFrame of coding alignments produced by BLAT.

    consensus metric w/rnaseq = 0.05 * cov + 0.5 * identity + 0.3 * rnaseq_support + 0.15 * structure score
    consensus metric w/o rnaseq = 0.05 * cov + 0.7 * identity + 0.25 * structure score

    structure score: how much of the parental transcript structure did we reproduce? Treat single exon transcripts
    as having no missing introns.
    structure score = 0.7 * percent original introns + 0.2 * percent original exons + 0.1 * evaluation score

    evaluation score: provides a bonus to transcripts without frameshifts and with proper CDS ends
    evaluation score = 1 - (I(in frame stop)
                            + I(coding indel)
                            + I(CdsStartStat = imcpl and == cmpl in reference)
                            + I(CdsEndStat = imcpl and == cmpl in reference)) / 4

    TranscriptClass: Classes were defined based on distribution fits previously as Passing/Failing. We extend that now
    to include the concept of Excellent, which are transcripts who were Passing previously and are now either:
    w/rnaseq: evaluation score == 1 and rnaseq_support > 0.8
    w/o rnaseq: evaluation score == 1
    """
    def evaluation_score(s):
        ifs = s.InFrameStop == 1
        indel = s.CodingDeletion == 0 or s.CodingInsertion == 0
        start_penalty = 1 if s.StartCodon_Tgt == 0 and s.StartCodon_Ref == 1 else 0
        stop_penalty = 1 if s.StopCodon_Tgt == 0 and s.StopCodon_Ref == 1 else 0
        # using bool -> integer implicit casting here
        return 1 - tools.mathOps.format_ratio(ifs + indel + start_penalty + stop_penalty, 4)

    def structure_score(s):
        # need to handle the case of a single exon transcript
        percent_original_introns = s.PercentOriginalIntrons if not math.isnan(s.PercentOriginalIntrons) else 1
        return 0.7 * percent_original_introns + 0.2 * percent_original_introns + 0.1 * evaluation_score(s)

    def rnaseq_support(s):
        if s.AlnMode == 'mRNA':
            r = s.PercentIntronsSupported
        else:
            r = s.PercentCodingIntronsSupported
        return 1 if np.isnan(r) else r

    def coding_consensus_metric(s):
        if has_rnaseq_data is True:
            return 0.05 * s.AlnCoverage + 0.5 * s.AlnIdentity + 0.3 * rnaseq_support(s) + 0.15 * structure_score(s)
        else:
            return 0.05 * s.AlnCoverage + 0.7 * s.AlnIdentity + 0.25 * structure_score(s)

    def calculate_class(s):
        """upgrades the class of a transcript if it was passing and now should be excellent"""
        if has_rnaseq_data is True:
            if s.TranscriptClass == 'Passing' and evaluation_score(s) == 1 and s.PercentIntronsSupported >= 0.8:
                return 'Excellent'
        else:
            if s.TranscriptClass == 'Passing' and evaluation_score(s) == 1:
                return 'Excellent'
        # TranscriptClass may start out as nan if this is a TM/TMR/CGP transcript
        if s.TranscriptClass == 'Failing' or s.TranscriptClass == 'Passing':
            return s.TranscriptClass
        return 'Failing'

    merged_df['ConsensusScore'] = merged_df.apply(coding_consensus_metric, axis=1)
    merged_df['TranscriptClass'] = merged_df.apply(calculate_class, axis=1)
    return merged_df


def score_non_coding_df(tm_eval, intron_df):
    """
    Scores a transMap DataFrame of alignments. Coding alignments will be removed.
    non coding consensus scores rely much more heavily on coverage, because these are only used in the case where
    all transcripts are failing, and we want the longest transcript of reasonable quality in that case.

    consensus metric w/rnaseq = 0.05 * identity + 0.5 * cov + 0.3 * rnaseq_support + 0.15 * percent original introns
    consensus metric w/o rnaseq = 0.05 * identity + 0.7 * cov + 0.25 * percent original introns

    TranscriptClass: A non-coding transcript can be upgraded if it has 80% intron support, which means upgrades can only
    occur if there was input RNA-seq data.

    :param tm_eval: DataFrame produced by load_transmap_evals()
    :param intron_df: DataFrame produced by load_intron_vectors(), if we had RNAseq data
    :return: DataFrame
    """
    def non_coding_consensus_metric(s):
        if intron_df is not None:
            return 0.05 * s.TransMapIdentity + 0.5 * s.TransMapCoverage + 0.15 * s.TransMapPercentOriginalIntrons + \
                   0.3 * s.PercentIntronsSupported if not math.isnan(s.PercentIntronsSupported) else 1
        else:
            return 0.05 * s.TransMapIdentity + 0.7 * s.TransMapCoverage + 0.25 * s.TransMapPercentOriginalIntrons

    def calculate_class(s):
        """upgrades the class of a transcript if it was passing and now should be excellent"""
        if s.TranscriptClass == 'Passing' and s.PercentIntronsSupported >= 0.8:
            return 'Excellent'
        return s.TranscriptClass

    df = tm_eval[tm_eval.TranscriptBiotype != 'protein_coding']
    # merge in intron support data, if we have it
    if intron_df is not None:
        # left outer join to keep transcripts not in intron dict (single exon)
        df = pd.merge(df, intron_df, on=['GeneId', 'AlignmentId'], how='left')
        df['TranscriptClass'] = df.apply(calculate_class, axis=1)

    df['ConsensusScore'] = df.apply(non_coding_consensus_metric, axis=1)
    # rename the identity/coverage columns to match the coding
    df = df.rename(columns={'TransMapIdentity': 'AlnIdentity', 'TransMapCoverage': 'AlnCoverage',
                            'TransMapBadness': 'AlnBadness'})
    df['AlnMode'] = ['mRNA'] * len(df)
    return df


def find_novel_transcripts(intron_df, tx_dict, metrics, num_introns=3, splice_support=0.8):
    """Finds novel transcripts, builds their attributes"""
    novel_transcripts = {tx.name: tx.name for tx in tx_dict.itervalues() if 'jg' in tx.name2}
    novel_df = intron_df[intron_df.AlignmentId.isin(novel_transcripts)]
    novel_df = novel_df[(novel_df.PercentIntronsSupported >= splice_support) & (novel_df.NumIntrons > num_introns)]
    novel_genes = set(novel_df.GeneId)
    metrics['CGP'] = {'Novel genes': len(novel_genes), 'Novel transcripts': len(novel_transcripts)}
    metrics['Transcript Modes']['augCGP'] += len(novel_transcripts)
    consensus = {}
    for novel_tx in novel_df.AlignmentId:
        consensus[novel_tx] = {'transcript_class': 'Novel', 'gene_biotype': 'unknown_likely_coding',
                               'transcript_biotype': 'unknown_likely_coding'}
    return consensus


def find_best_score(tx_df, column='ConsensusScore'):
    """
    Finds the best transcript in the pre-sorted filtered DataFrame, handling the case where it does not exist.
    """
    try:
        if isinstance(tx_df, pd.core.series.Series):
            return pd.DataFrame([tx_df])
        else:
            best_score = tx_df.iloc[0][column]
    except IndexError:
        return None
    return tx_df[tx_df[column] == best_score]


def incorporate_tx(best_rows, gene_id, metrics, failed_gene):
    """incorporate a transcript into the consensus set, storing metrics. Updates gene_seen."""
    best_series = best_rows.iloc[0]
    # construct the tags for this transcript
    d = {'source_transcript': best_series.name,
         'source_gene': gene_id,
         'transcript_mode': tools.nameConversions.alignment_type(best_series.AlignmentId),
         'score': round(best_series.ConsensusScore, 2),
         'failed_gene': failed_gene,
         'transcript_class': best_series.TranscriptClass,
         'gene_biotype': best_series.GeneBiotype,
         'transcript_biotype': best_series.TranscriptBiotype}
    if best_series.Paralogy > 1:
        assert best_series.ParalogStatus is not None
        d['paralogy'] = best_series.Paralogy
        d['paralog_status'] = best_series.ParalogStatus
    if best_series.GeneAlternateContigs is not None:
        d['gene_alterate_contigs'] = best_series.GeneAlternateContigs
    if best_series.GeneName is not None:
        d['source_gene_common_name'] = best_series.GeneName
    if best_series.TranscriptBiotype == 'protein_coding':
        metrics['Transcript Modes'][evaluate_ties(best_rows)] += 1
    metrics['Transcript Categories'][best_series.TranscriptBiotype][best_series.TranscriptClass] += 1
    metrics['Coverage'][best_series.TranscriptBiotype].append(best_series.AlnCoverage)
    metrics['Identity'][best_series.TranscriptBiotype].append(best_series.AlnIdentity)
    metrics['Consensus Score'][best_series.TranscriptBiotype].append(best_series.ConsensusScore)
    metrics['Splice Support'][best_series.TranscriptBiotype].append(best_series.PercentIntronsSupported)
    return best_series.AlignmentId, d


def evaluate_ties(best_rows):
    """Find out how many transcript modes agreed on this"""
    return ','.join(sorted(set([tools.nameConversions.alignment_type(x) for x in best_rows.AlignmentId])))


def is_failed_df(df):
    """Failed genes have no passing/excellent transcripts. Handles series"""
    try:
        return not (df.TranscriptClass.str.contains('Passing').any() or
                    df.TranscriptClass.str.contains('Excellent').any())
    except AttributeError:
        return df.TranscriptClass == 'Failing'


def rescue_failed_gene(gene_df, tx_dict, gene_id, metrics):
    """Rescues a failed gene by picking the one transcript with the highest length * score"""
    tx_lengths = np.array([len(tx_dict[x]) for x in gene_df.AlignmentId])
    gene_df['RescueScore'] = gene_df.ConsensusScore * tx_lengths
    gene_df = gene_df.sort_values('RescueScore', ascending=False)
    best_rows = find_best_score(gene_df, 'RescueScore')
    return incorporate_tx(best_rows, gene_id, metrics, failed_gene=False)


def slice_df(df, ix):
    """
    Slices a DataFrame by an index, handling the case where the index is missing
    """
    try:
        return df.xs(ix)
    except KeyError:
        return pd.DataFrame()


def find_novel_cgp_splices(gene_consensus_dict, gene_df, tx_dict, gene_id, common_name_map, metrics, failed_gene):
    """
    Finds novel splice junctions in CGP transcripts. If there are any, these get included as a novel isoform.
    """
    existing_splices = set()
    for consensus_tx in gene_consensus_dict:
        existing_splices.update(tx_dict[consensus_tx].intron_intervals)

    # extract CGP transcripts and slice the intron vector
    # this will collapse the duplicate vectors into one
    cgp_df = gene_df[tools.nameConversions.aln_id_is_cgp(gene_df.AlignmentId.str)]
    cgp_introns = {s.AlignmentId: s.IntronVector for _, s in cgp_df.iterrows()}
    cgp_tx_dict = {}
    for cgp_tx, intron_vector in cgp_introns.iteritems():
        cgp_tx_obj = tx_dict[cgp_tx]
        for interval, intron_score in zip(*[cgp_tx_obj.intron_intervals, intron_vector]):
            if intron_score > 0 and interval not in existing_splices:
                metrics['Novel isoforms'] += 1
                metrics['Transcript Modes']['augCGP'] += 1
                cgp_tx_dict[cgp_tx] = {'transcript_class': 'novel', 'source_gene': gene_id,
                                       'failed_gene': failed_gene, 'transcript_mode': 'augCGP',
                                       'transcript_biotype': 'unknown_likely_coding',
                                       'gene_biotype': 'unknown_likely_coding'}
                common_name = common_name_map[gene_id]
                if common_name != gene_id:
                    cgp_tx_dict[cgp_tx]['source_gene_common_name'] = common_name
    return cgp_tx_dict


def deduplicate_consensus(consensus_dict, tx_dict, metrics):
    """
    In the process of consensus building, we may find that we have ended up with more than one transcript for a gene
    that are actually identical. Remove these, picking the best based on their score, favoring the transcript
    whose biotype matches the parent.
    """
    def resolve_duplicate(tx_list, consensus_dict):
        biotype_txs = [tx for tx in tx_list if
                       consensus_dict[tx].get('gene_biotype', None) == consensus_dict[tx].get('transcript_biotype', None)]
        if len(biotype_txs) > 0:
            sorted_scores = sorted([[tx, consensus_dict[tx].get('score', 0)]for tx in biotype_txs], key=lambda (tx, s): -s)
            return sorted_scores[0][0]
        else:
            sorted_scores = sorted([[tx, consensus_dict[tx].get('score', 0)]for tx in tx_list], key=lambda (tx, s): -s)
            return sorted_scores[0][0]

    def add_duplicate_field(best_tx, tx_list, consensus_dict, deduplicated_consensus):
        deduplicated_consensus[best_tx] = consensus_dict[best_tx]
        tx_list = [tools.nameConversions.strip_alignment_numbers(aln_id) for aln_id in tx_list]
        best_tx_base = tools.nameConversions.strip_alignment_numbers(best_tx)
        deduplicated_consensus[best_tx]['alternative_source_transcripts'] = ','.join(set(tx_list) - {best_tx_base})

    # build a dictionary mapping duplicates making use of hashing intervals
    duplicates = collections.defaultdict(list)
    for aln_id in consensus_dict:
        tx = tx_dict[aln_id]
        duplicates[frozenset(tx.exon_intervals)].append(aln_id)

    # begin iterating
    deduplicated_consensus = {}
    for tx_list in duplicates.itervalues():
        if len(tx_list) > 1:
            metrics['Duplicate transcripts'][len(tx_list)] += 1
            best_tx = resolve_duplicate(tx_list, consensus_dict)
            add_duplicate_field(best_tx, tx_list, consensus_dict, deduplicated_consensus)
        else:
            tx_id = tx_list[0]
            deduplicated_consensus[tx_id] = consensus_dict[tx_id]

    return deduplicated_consensus


def resolve_opposite_strand(deduplicated_consensus, tx_dict, metrics):
    """
    Resolves situations where multiple transcripts of the same gene are on opposite strands. Does so by looking for
    the largest sum of scores.
    """
    gene_dict = collections.defaultdict(list)
    for tx_id, attrs in deduplicated_consensus.iteritems():
        tx_obj = tx_dict[tx_id]
        gene_dict[tx_obj.name2].append([tx_obj, attrs])

    deduplicated_strand_resolved_consensus = []
    for gene in gene_dict:
        tx_objs, attrs = zip(*gene_dict[gene])
        if len(set(tx_obj.strand for tx_obj in tx_objs)) > 1:
            strand_scores = collections.Counter()
            for tx_obj, attrs in gene_dict[gene]:
                strand_scores[tx_obj.strand] += attrs.get('score', 0)
            best_strand = sorted(strand_scores.items())[0][0]
            for tx_obj, attrs in gene_dict[gene]:
                if tx_obj.strand == best_strand:
                    deduplicated_strand_resolved_consensus.append([tx_obj.name, attrs])
                else:
                    metrics['Discarded by strand resolution'] += 1
        else:
            deduplicated_strand_resolved_consensus.extend([[tx_obj.name, attrs] for tx_obj, attrs in gene_dict[gene]])
    return deduplicated_strand_resolved_consensus


def calculate_completeness(final_consensus, metrics):
    """calculates final completeness to make arithmetic easier"""
    genes = collections.defaultdict(set)
    txs = collections.Counter()
    for aln_id, c in final_consensus:
        if c['transcript_biotype'] == 'unknown_likely_coding':
            continue
        genes[c['gene_biotype']].add(c['source_gene'])
        txs[c['transcript_biotype']] += 1
    genes = {biotype: len(gene_list) for biotype, gene_list in genes.iteritems()}
    metrics['Completeness'] = {'Gene': genes, 'Transcript': txs}


def write_consensus_gps(consensus_gp, consensus_gp_info, final_consensus, tx_dict, genome):
    """
    Write the resulting gp + gp_info, generating genome-specific unique identifiers
    """
    gene_count = 0
    tx_count = 1
    consensus_gene_dict = DefaultOrderedDict(lambda: DefaultOrderedDict(list))  # used to make gff3 next
    gp_infos = []
    genes_seen = set()
    consensus_gp = luigi.LocalTarget(consensus_gp)
    with consensus_gp.open('w') as out_gp:
        for tx, attrs in final_consensus:
            tx_obj = copy.deepcopy(tx_dict[tx])
            tx_obj.name = id_template.format(genome=genome, tag_type='T', unique_id=tx_count)
            tx_obj.id = attrs.get('score', 0)
            tx_count += 1
            if tx_obj.name2 not in genes_seen:
                genes_seen.add(tx_obj.name2)
                gene_count += 1
            tx_obj.name2 = id_template.format(genome=genome, tag_type='G', unique_id=gene_count)
            out_gp.write('\t'.join(tx_obj.get_gene_pred()) + '\n')
            consensus_gene_dict[tx_obj.chromosome][tx_obj.name2].append([tx_obj, attrs.copy()])
            gp_info = attrs.copy()
            gp_info['transcript_id'] = tx_obj.name
            gp_info['gene_id'] = tx_obj.name2
            gp_infos.append(gp_info)
    gp_info_df = pd.DataFrame(gp_infos)
    gp_info_df = gp_info_df.set_index(['gene_id', 'transcript_id'])
    consensus_gp_info = luigi.LocalTarget(consensus_gp_info)
    with consensus_gp_info.open('w') as outf:
        gp_info_df.to_csv(outf, sep='\t')
    return consensus_gene_dict


def write_consensus_gff3(consensus_gene_dict, consensus_gff3):
    """
    Write the consensus set in gff3 format
    """
    def convert_frame(exon_frame):
        """converts genePred-style exonFrame to GFF-style phase"""
        mapping = {0: 0, 1: 2, 2: 1, -1: '.'}
        return mapping[exon_frame]

    def convert_attrs(attrs, id_field):
        """converts the attrs dict to a attributes field. assigns name to the gene common name for display"""
        attrs['ID'] = id_field
        try:
            score = attrs['score']
            del attrs['score']
        except KeyError:
            score = 0
        if 'source_gene_common_name' in attrs:
            attrs['Name'] = attrs['source_gene_common_name']
        attrs_str = ['='.join([key, str(val)]) for key, val in sorted(attrs.iteritems())]
        return score, ';'.join(attrs_str)

    def generate_gene_record(chrom, tx_objs, gene_id, attrs):
        """calculates the gene interval for this list of tx"""
        intervals = set()
        for tx in tx_objs:
            intervals.update(tx.exon_intervals)
        intervals = sorted(intervals)
        strand = tx_objs[0].strand
        # subset the attrs to gene fields
        useful_keys = ['source_gene_common_name', 'source_gene', 'gene_biotype', 'failed_gene',
                       'alternative_source_transcripts', 'paralog_status', 'gene_alterate_contigs']
        attrs = {key: attrs[key] for key in useful_keys if key in attrs}
        score, attrs_field = convert_attrs(attrs, gene_id)
        return [chrom, 'CAT', 'gene', intervals[0].start + 1, intervals[-1].stop + 1, score, strand, '.', attrs_field]

    def generate_transcript_record(chrom, tx_obj, attrs):
        """generates transcript records, calls generate_exon_records to generate those too"""
        tx_id = tx_obj.name
        gene_id = tx_obj.name2
        attrs['Parent'] = gene_id
        score, attrs_field = convert_attrs(attrs, tx_id)
        yield [chrom, 'CAT', 'transcript', tx_obj.start + 1, tx_obj.stop + 1, score, tx_obj.strand, '.', attrs_field]
        for line in generate_exon_records(chrom, tx_obj, tx_id, attrs):
            yield line
        for line in generate_start_stop_codon_records(chrom, tx_obj, tx_id, attrs):
            yield line

    def generate_exon_records(chrom, tx_obj, tx_id, attrs):
        """generates exon records"""
        attrs['Parent'] = tx_id
        for i, (exon, exon_frame) in enumerate(zip(*[tx_obj.exon_intervals, tx_obj.exon_frames]), 1):
            score, attrs_field = convert_attrs(attrs, 'exon:{}:{}'.format(tx_id, i))
            yield [chrom, 'CAT', 'exon', exon.start + 1, exon.stop + 1, score, exon.strand, '.', attrs_field]
            cds_interval = exon.intersection(tx_obj.coding_interval)
            if cds_interval is not None:
                score, attrs_field = convert_attrs(attrs, 'CDS:{}:{}'.format(tx_id, i))
                yield [chrom, 'CAT', 'CDS', cds_interval.start + 1, cds_interval.stop + 1, score, exon.strand,
                       convert_frame(exon_frame), attrs_field]

    def generate_start_stop_codon_records(chrom, tx_obj, tx_id, attrs):
        """generate start/stop codon GFF3 records, handling frame appropriately"""
        cds_frames = [x for x in tx_obj.exon_frames if x != -1]
        if tx_obj.cds_start_stat == 'cmpl':
            score, attrs_field = convert_attrs(attrs, 'start_codon:{}'.format(tx_id))
            start, stop = tools.transcripts.get_start_interval(tx_obj)
            if tx_obj.strand == '-':
                start_frame = convert_frame(cds_frames[-1])
            else:
                start_frame = convert_frame(cds_frames[0])
            yield [chrom, 'CAT', 'start_codon', start + 1, stop + 1, score, tx_obj.strand, start_frame, attrs_field]
        if tx_obj.cds_end_stat == 'cmpl':
            score, attrs_field = convert_attrs(attrs, 'stop_codon:{}'.format(tx_id))
            start, stop = tools.transcripts.get_stop_interval(tx_obj)
            if tx_obj.strand == '-':
                stop_frame = convert_frame(cds_frames[-1])
            else:
                stop_frame = convert_frame(cds_frames[0])
            yield [chrom, 'CAT', 'stop_codon', start + 1, stop + 1, score, tx_obj.strand, stop_frame, attrs_field]

    # main gff3 writing logic
    consensus_gff3 = luigi.LocalTarget(consensus_gff3)
    with consensus_gff3.open('w') as out_gff3:
        out_gff3.write('##gff-version 3\n')
        for chrom in sorted(consensus_gene_dict):
            for gene_id, tx_list in consensus_gene_dict[chrom].iteritems():
                tx_objs, attrs_list = zip(*tx_list)
                attrs = tx_list[0][1]  # grab the attrs from the first transcript
                tools.fileOps.print_row(out_gff3, generate_gene_record(chrom, tx_objs, gene_id, attrs))
                tx_lines = []
                for tx_obj, attrs in tx_list:
                    tx_lines.extend(list(generate_transcript_record(chrom, tx_obj, attrs)))
                tx_lines = sorted(tx_lines, key=lambda l: l[3])
                tools.fileOps.print_rows(out_gff3, tx_lines)
