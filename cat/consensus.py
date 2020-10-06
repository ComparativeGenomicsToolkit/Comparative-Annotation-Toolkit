"""
Generates consensus gene set.

This module takes as input the genePreds produced by transMap, AugustusTM(R) and AugustusCGP and generates a consensus
of these, producing a filtered gene set.

This process relies on a combination of metrics and evaluations loaded to a sqlite database by the classify module.


GFF3 tags generated in this process:
1. source_transcript: The name of the parent transcript, if it exists
2. source_gene: The name of the parent gene, if it exists
3. source_gene_common_name: The common name of the parent gene, if it is different from the source gene
4. transcript_mode: The name of the mode of operation that generated this transcript
5. transcript_class: One of possible_paralog, poor_alignment, putative_novel, putative_novel_isoform, ortholog
6. paralogy: The names of paralogous alignments
7. unfiltered_paralogy: The names of all alignments without globalNearBest filtering.
8. gene_biotype: gene biotype
9. transcript_biotype: transcript biotype
10. alternative_source_transcripts: Other possible transcripts, if this was collapsed as the result of deduplication
11. gene_alternate_contigs: contigs that this gene was also found on are comma separated in this tag.
12: transcript_modes: The mode(s) that generated this transcript
"""
import collections
import luigi
import logging
import pandas as pd

import tools.bio
import tools.intervals
import tools.misc
import tools.mathOps
import tools.fileOps
import tools.sqlInterface
import tools.transcripts
import tools.nameConversions
import tools.procOps
from tools.defaultOrderedDict import DefaultOrderedDict

logger = logging.getLogger('cat')

id_template = '{genome:.10}_{tag_type}{unique_id:07d}'


def generate_consensus(args):
    """
    Main consensus finding logic.

    :param args: Argument namespace from luigi
    """
    # load all genePreds
    tx_dict = tools.transcripts.load_gps(args.gp_list)
    # load reference annotation information
    ref_df = tools.sqlInterface.load_annotation(args.ref_db_path)
    ref_biotype_counts = collections.Counter(ref_df.TranscriptBiotype)
    coding_count = ref_biotype_counts['protein_coding']
    non_coding_count = sum(y for x, y in ref_biotype_counts.items() if x != 'protein_coding')
    # gene transcript map to iterate over so that we capture missing gene information
    gene_biotype_map = tools.sqlInterface.get_gene_biotype_map(args.ref_db_path)
    transcript_biotype_map = tools.sqlInterface.get_transcript_biotype_map(args.ref_db_path)
    # load transMap evaluation data
    tm_eval_df = load_transmap_evals(args.db_path)
    # load the homGeneMapping data for transMap/augTM/augTMR
    tx_modes = [x for x in args.tx_modes if x in ['transMap', 'augTM', 'augTMR']]
    hgm_df = pd.concat([load_hgm_vectors(args.db_path, tx_mode) for tx_mode in tx_modes])
    # load the alignment metrics data
    mrna_metrics_df = pd.concat([load_metrics_from_db(args.db_path, tx_mode, 'mRNA') for tx_mode in tx_modes])
    cds_metrics_df = pd.concat([load_metrics_from_db(args.db_path, tx_mode, 'CDS') for tx_mode in tx_modes])
    eval_df = pd.concat([load_evaluations_from_db(args.db_path, tx_mode) for tx_mode in tx_modes]).reset_index()
    coding_df, non_coding_df = combine_and_filter_dfs(tx_dict, hgm_df, mrna_metrics_df, cds_metrics_df, tm_eval_df,
                                                      ref_df, eval_df, args.intron_rnaseq_support,
                                                      args.exon_rnaseq_support, args.intron_annot_support,
                                                      args.exon_annot_support, args.original_intron_support,
                                                      args.in_species_rna_support_only)
    if len(coding_df) + len(non_coding_df) == 0:
        raise RuntimeError('No transcripts pass filtering for species {}. '
                           'Consider lowering requirements. Please see the manual.'.format(args.genome))
    elif len(coding_df) == 0 and coding_count > 0:
        logger.warning('No protein coding transcripts pass filtering for species {}. '
                       'Consider lowering requirements. Please see the manual.'.format(args.genome))
    elif len(non_coding_df) == 0 and non_coding_count > 0:
        logger.warning('No non-coding transcripts pass filtering for species {}. '
                       'Consider lowering requirements. Please see the manual.'.format(args.genome))
    scored_coding_df, scored_non_coding_df = score_filtered_dfs(coding_df, non_coding_df,
                                                                args.in_species_rna_support_only)
    scored_df = merge_scored_dfs(scored_coding_df, scored_non_coding_df)
    best_alignments = scored_df.groupby('TranscriptId')['TranscriptScore'].transform(max) == scored_df['TranscriptScore']
    best_df = scored_df[best_alignments].reset_index()

    # store some metrics for plotting
    metrics = {'Transcript Missing': collections.Counter(),
               'Gene Missing': collections.Counter(),
               'Transcript Modes': collections.Counter(),  # coding only
               'Duplicate transcripts': collections.Counter(),
               'Discarded by strand resolution': 0,
               'Coverage': collections.defaultdict(list),
               'Identity': collections.defaultdict(list),
               'Splice Support': collections.defaultdict(list),
               'Exon Support': collections.defaultdict(list),
               'Original Introns': collections.defaultdict(list),
               'Splice Annotation Support': collections.defaultdict(list),
               'Exon Annotation Support': collections.defaultdict(list),
               'IsoSeq Transcript Validation': collections.Counter()}

    # we can keep track of missing stuff now
    for gene_biotype, tx_df in best_df.groupby('GeneBiotype'):
        biotype_genes = {gene_id for gene_id, b in gene_biotype_map.items() if b == gene_biotype}
        metrics['Gene Missing'][gene_biotype] = len(biotype_genes) - len(set(tx_df.GeneId))
    for tx_biotype, tx_df in best_df.groupby('TranscriptBiotype'):
        biotype_txs = {gene_id for gene_id, b in transcript_biotype_map.items() if b == tx_biotype}
        metrics['Transcript Missing'][tx_biotype] = len(biotype_txs) - len(set(tx_df.TranscriptId))

    # main consensus finding -- using incorporate_tx to transform best scoring transcripts
    # stores a mapping of alignment IDs to tags for the final consensus set
    consensus_dict = {}
    for (gene_id, tx_id), s in best_df.groupby(['GeneId', 'TranscriptId']):
        aln_id, m = incorporate_tx(s, gene_id, metrics, args.hints_db_has_rnaseq)
        consensus_dict[aln_id] = m

    # if we ran in either denovo mode, load those data and detect novel genes
    if len(args.denovo_tx_modes) > 0:
        metrics['denovo'] = {}
        for tx_mode in args.denovo_tx_modes:
            metrics['denovo'][tx_mode] = {'Possible paralog': 0, 'Poor alignment': 0, 'Putative novel': 0,
                                          'Possible fusion': 0, 'Putative novel isoform': 0, 'User filtered': 0}
        denovo_dict = find_novel(args.db_path, tx_dict, consensus_dict, ref_df, metrics, gene_biotype_map,
                                 args.denovo_num_introns, args.in_species_rna_support_only,
                                 args.denovo_tx_modes, args.denovo_splice_support, args.denovo_exon_support,
                                 args.denovo_ignore_novel_genes, args.denovo_novel_end_distance,
                                 args.denovo_allow_unsupported, args.denovo_allow_bad_annot_or_tm,
                                 args.denovo_only_novel_genes, args.denovo_allow_novel_ends)
        consensus_dict.update(denovo_dict)

    # perform final filtering steps
    deduplicated_consensus = deduplicate_consensus(consensus_dict, tx_dict, metrics)
    deduplicated_strand_resolved_consensus = resolve_opposite_strand(deduplicated_consensus, tx_dict, metrics)

    if 'augPB' in args.denovo_tx_modes:
        deduplicated_strand_resolved_consensus = validate_pacbio_splices(deduplicated_strand_resolved_consensus,
                                                                         args.db_path, tx_dict, metrics,
                                                                         args.require_pacbio_support)

    if args.filter_overlapping_genes is True:
        gene_resolved_consensus = resolve_overlapping_cds_intervals(args.overlapping_gene_distance,
                                                                    deduplicated_strand_resolved_consensus, tx_dict)
    else:
        gene_resolved_consensus = deduplicated_strand_resolved_consensus

    # sort by genomic interval for prettily increasing numbers
    final_consensus = sorted(gene_resolved_consensus,
                             key=lambda tx_attrs: (tx_dict[tx_attrs[0]].chromosome, tx_dict[tx_attrs[0]].start))

    # calculate final gene set completeness
    calculate_completeness(final_consensus, metrics)
    # add some interesting metrics on how much using Augustus modes improved our results
    if 'augTM' in tx_modes or 'augTMR' in tx_modes:
        calculate_improvement_metrics(final_consensus, scored_df, tm_eval_df, hgm_df, metrics)
    calculate_indel_metrics(final_consensus, eval_df, metrics)
    # write out results. consensus tx dict has the unique names
    consensus_gene_dict = write_consensus_gps(args.consensus_gp, args.consensus_gp_info,
                                              final_consensus, tx_dict, args.genome)
    write_consensus_gff3(consensus_gene_dict, args.consensus_gff3)
    write_consensus_fastas(consensus_gene_dict, args.consensus_fasta, args.consensus_protein_fasta, args.fasta)

    return metrics


def load_transmap_evals(db_path):
    """
    Loads the database tables associated with both transMap filtering and transMap evaluation, merging them. Keep only
    the columns related to paralog resolution and classification.
    """
    # load transMap results
    tm_eval = tools.sqlInterface.load_alignment_evaluation(db_path)
    # load transMap filtering results
    tm_filter_eval = tools.sqlInterface.load_filter_evaluation(db_path)

    # combine transMap evaluation and transMap filtering into one table
    # the transMap filtering columns are used for tags in the output
    tm_eval_df = pd.merge(tm_eval, tm_filter_eval, on=['TranscriptId', 'AlignmentId'])
    return tm_eval_df.drop('AlignmentId', axis=1)


def calculate_vector_support(s, resolve_nan=None, num_digits=4):
    """For vectors parsed by parse_text_vector(), convert to a percentage between 0 and 100"""
    return 100 * tools.mathOps.format_ratio(len([x for x in s if x > 0]), len(s), resolve_nan=resolve_nan,
                                            num_digits=num_digits)


def load_hgm_vectors(db_path, tx_mode):
    """
    Loads the intron vector table output by the homGeneMapping module. Returns a DataFrame with the parsed vectors
    as well as the combined score based on tx_mode -- for augCGP we look at the CDS annotation score rather than the
    exon score because CGP has a coding-only model.
    """
    session = tools.sqlInterface.start_session(db_path)
    intron_table = tools.sqlInterface.tables['hgm'][tx_mode]
    hgm_df = tools.sqlInterface.load_intron_vector(intron_table, session)

    # start calculating support levels for consensus finding
    cols = ['IntronAnnotSupport', 'ExonAnnotSupport', 'CdsAnnotSupport',
            'ExonRnaSupport', 'IntronRnaSupport',
            'AllSpeciesExonRnaSupport', 'AllSpeciesIntronRnaSupport']
    for col in cols:
        hgm_df[col] = [list(map(int, x)) if len(x[0]) > 0 else [] for x in hgm_df[col].str.split(',').tolist()]
        hgm_df[col + 'Percent'] = hgm_df[col].apply(calculate_vector_support, resolve_nan=1)
    return hgm_df


def load_metrics_from_db(db_path, tx_mode, aln_mode):
    """
    Loads the alignment metrics for the mRNA/CDS alignments of transMap/AugustusTM/TMR
    """
    session = tools.sqlInterface.start_session(db_path)
    metrics_table = tools.sqlInterface.tables[aln_mode][tx_mode]['metrics']
    metrics_df = tools.sqlInterface.load_metrics(metrics_table, session)
    # unstack flattens the long-form data structure
    metrics_df = metrics_df.set_index(['AlignmentId', 'classifier']).unstack('classifier')
    metrics_df.columns = [col[1] for col in metrics_df.columns]
    metrics_df = metrics_df.reset_index()
    cols = ['AlnCoverage', 'AlnGoodness', 'AlnIdentity', 'PercentUnknownBases']
    metrics_df[cols] = metrics_df[cols].apply(pd.to_numeric)
    metrics_df['OriginalIntrons'] = metrics_df['OriginalIntrons'].fillna('')
    metrics_df['OriginalIntrons'] = [list(map(int, x)) if len(x[0]) > 0 else [] for x in
                                     metrics_df['OriginalIntrons'].str.split(',').tolist()]
    metrics_df['OriginalIntronsPercent'] = metrics_df['OriginalIntrons'].apply(calculate_vector_support, resolve_nan=1)
    session.close()
    return metrics_df


def load_evaluations_from_db(db_path, tx_mode):
    """
    Loads the indel information from the evaluation database. We give preference to CDS alignments, but fall back
    to mRNA alignments.
    """
    def aggfunc(s):
        """
        Preferentially pick CDS stats over mRNA stats, if they exist
        They only exist for coding transcripts, and only those whose CDS alignments didn't fail
        """
        if s.value_CDS.any():
            c = set(s[s.value_CDS > 0].name)
        else:
            c = set(s[s.value_mRNA > 0].name)
        cols = ['Frameshift', 'CodingInsertion', 'CodingDeletion', 'CodingMult3Indel']
        return pd.Series(('CodingDeletion' in c or 'CodingInsertion' in c,
                          'CodingInsertion' in c, 'CodingDeletion' in c,
                          'CodingMult3Deletion' in c or 'CodingMult3Insertion' in c), index=cols)

    session = tools.sqlInterface.start_session(db_path)
    cds_table = tools.sqlInterface.tables['CDS'][tx_mode]['evaluation']
    mrna_table = tools.sqlInterface.tables['mRNA'][tx_mode]['evaluation']
    cds_df = tools.sqlInterface.load_evaluation(cds_table, session)
    mrna_df = tools.sqlInterface.load_evaluation(mrna_table, session)
    cds_df = cds_df.set_index('AlignmentId')
    mrna_df = mrna_df.set_index('AlignmentId')
    merged = mrna_df.reset_index().merge(cds_df.reset_index(), how='outer', on=['AlignmentId', 'name'],
                                         suffixes=['_mRNA', '_CDS'])
    eval_df = merged.groupby('AlignmentId').apply(aggfunc)
    return eval_df


def load_alt_names(db_path, denovo_tx_modes):
    """Load the alternative tx tables for augCGP/augPB"""
    session = tools.sqlInterface.start_session(db_path)
    r = []
    for tx_mode in denovo_tx_modes:
        table = tools.sqlInterface.tables['alt_names'][tx_mode]
        r.append(tools.sqlInterface.load_alternatives(table, session))
    df = pd.concat(r)
    # rename TranscriptId to AlignmentId. This is all super confusing and silly
    # the reason is that homGeneMapping has the gene -> tx -> aln ID hierarchy we inherit to simplify things
    df.columns = [x if x != 'TranscriptId' else 'AlignmentId' for x in df.columns]
    return df


def combine_and_filter_dfs(tx_dict, hgm_df, mrna_metrics_df, cds_metrics_df, tm_eval_df, ref_df, eval_df,
                           intron_rnaseq_support, exon_rnaseq_support, intron_annot_support, exon_annot_support,
                           original_intron_support, in_species_rna_support_only):
    """
    Updates the DataFrame based on support levels. Filters based on user-tunable flags for support levels.
    :param tx_dict: dictionary of genePredTranscript objects. Used to remove things filtered out by transMap
    :param hgm_df: df produced by load_hgm_vectors() (all transcripts)
    :param mrna_metrics_df: df produced by load_metrics() (coding transcripts only) for mRNA alignments
    :param cds_metrics_df: df produced by load_metrics() (coding transcripts only) for CDS alignments
    :param tm_eval_df: df produced by load_transmap_evals() (all transcripts)
    :param ref_df: df produced by tools.sqlInterface.load_annotation()
    :param eval_df: produced by load_evaluations_from_db() (coding transcripts only)
    :param intron_rnaseq_support: Value 0-100. Percent of introns that must be supported by RNAseq
    :param exon_rnaseq_support: Value 0-100. Percent of exons supported by RNA-seq.
    :param intron_annot_support: Value 0-100. Percent of introns supported by the reference.
    :param exon_annot_support: Value 0-100. Percent of exons supported by the reference.
    :param original_intron_support: Value 0-100. Percent of introns that must be supported by this specific annotation.
    :param in_species_rna_support_only: Should we use the homGeneMapping vectors within-species or all-species?
    :return: filtered and merged dataframe
    """
    # add the reference information to gain biotype information
    hgm_ref_df = pd.merge(hgm_df, ref_df, on=['GeneId', 'TranscriptId'])
    # combine in homGeneMapping results
    hgm_ref_tm_df = pd.merge(hgm_ref_df, tm_eval_df, on=['GeneId', 'TranscriptId'])
    # remove filtered transMap
    hgm_ref_tm_df = hgm_ref_tm_df[hgm_ref_tm_df.AlignmentId.isin(tx_dict.keys())]
    # split merged_df into coding and noncoding
    coding_df = hgm_ref_tm_df[hgm_ref_tm_df.TranscriptBiotype == 'protein_coding']
    non_coding_df = hgm_ref_tm_df[hgm_ref_tm_df.TranscriptBiotype != 'protein_coding']
    # add metrics information to coding df
    metrics_df = pd.merge(mrna_metrics_df, cds_metrics_df, on='AlignmentId', suffixes=['_mRNA', '_CDS'])
    coding_df = pd.merge(coding_df, metrics_df, on='AlignmentId')
    # add evaluation information to coding df, where possible. This adds information on frame shifts.
    coding_df = pd.merge(coding_df, eval_df, on='AlignmentId', how='left')
    # fill the original intron values to 100 so we don't filter them out -- means a no-intron gene
    coding_df['OriginalIntronsPercent_mRNA'] = coding_df.OriginalIntronsPercent_mRNA.fillna(100)
    coding_df['OriginalIntronsPercent_CDS'] = coding_df.OriginalIntronsPercent_CDS.fillna(100)
    non_coding_df['TransMapOriginalIntronsPercent'] = non_coding_df.TransMapOriginalIntronsPercent.fillna(100)

    # huge ugly filtering expression for coding transcripts
    if in_species_rna_support_only is True:
        filt = ((coding_df.OriginalIntronsPercent_mRNA >= original_intron_support) &
                (coding_df.IntronAnnotSupportPercent >= intron_annot_support) &
                (coding_df.IntronRnaSupportPercent >= intron_rnaseq_support) &
                (coding_df.ExonAnnotSupportPercent >= exon_annot_support) &
                (coding_df.ExonRnaSupportPercent >= exon_rnaseq_support))
    else:
        filt = ((coding_df.OriginalIntronsPercent_mRNA >= original_intron_support) &
                (coding_df.IntronAnnotSupportPercent >= intron_annot_support) &
                (coding_df.AllSpeciesIntronRnaSupportPercent >= intron_rnaseq_support) &
                (coding_df.ExonAnnotSupportPercent >= exon_annot_support) &
                (coding_df.AllSpeciesExonRnaSupportPercent >= exon_rnaseq_support))
    coding_df = coding_df[filt]

    # huge ugly filtering expression for non coding transcripts
    if in_species_rna_support_only is True:
        filt = ((non_coding_df.TransMapOriginalIntronsPercent >= original_intron_support) &
                (non_coding_df.IntronAnnotSupportPercent >= intron_annot_support) &
                (non_coding_df.IntronRnaSupportPercent >= intron_rnaseq_support) &
                (non_coding_df.ExonAnnotSupportPercent >= exon_annot_support) &
                (non_coding_df.ExonRnaSupportPercent >= exon_rnaseq_support))
    else:
        filt = ((non_coding_df.TransMapOriginalIntronsPercent >= original_intron_support) &
                (non_coding_df.IntronAnnotSupportPercent >= intron_annot_support) &
                (non_coding_df.AllSpeciesIntronRnaSupportPercent >= intron_rnaseq_support) &
                (non_coding_df.ExonAnnotSupportPercent >= exon_annot_support) &
                (non_coding_df.AllSpeciesExonRnaSupportPercent >= exon_rnaseq_support))
    non_coding_df = non_coding_df[filt]

    return coding_df, non_coding_df


def score_filtered_dfs(coding_df, non_coding_df, in_species_rna_support_only):
    """
    Scores the alignments. The score is the additive combination of the following features:
    1) Alignment identity.
    2) Alignment coverage.
    3) Intron annotation support.
    4) Exon annotation support.
    5) Original intron support.
    If we have RNA-seq data, the following fields are also incorporated:
    6) Intron RNA-seq support.
    7) Exon RNA-seq support.

    In some future world these features could be combined differently, maybe with some fancy machine learning.

    Returns the dataframe sorted by scores after indexing.
    """
    def score(s):
        aln_id = s.AlnIdentity_CDS if s.TranscriptBiotype == 'protein_coding' else s.TransMapIdentity
        aln_cov = s.AlnCoverage_CDS if s.TranscriptBiotype == 'protein_coding' else s.TransMapCoverage
        orig_intron = s.OriginalIntronsPercent_mRNA if s.TranscriptBiotype == 'protein_coding' else s.TransMapOriginalIntronsPercent
        if in_species_rna_support_only:
            rna_support = s.ExonRnaSupportPercent + s.IntronRnaSupportPercent
        else:
            rna_support = s.AllSpeciesExonRnaSupportPercent + s.AllSpeciesIntronRnaSupportPercent
        return aln_id + aln_cov + s.IntronAnnotSupportPercent + s.ExonAnnotSupportPercent + orig_intron + rna_support

    for df in [coding_df, non_coding_df]:
        if len(df) > 0:
            df['TranscriptScore'] = df.apply(score, axis=1)
    return coding_df, non_coding_df


def merge_scored_dfs(scored_coding_df, scored_non_coding_df):
    """Merges the scored dataframes by changing some names around"""
    # for every non-coding TransMap metric, copy it to the other name
    for m in ['Coverage', 'Identity', 'Goodness']:
        scored_non_coding_df['Aln' + m + '_mRNA'] = scored_non_coding_df['TransMap' + m]
    merged_df = pd.concat([scored_non_coding_df, scored_coding_df])
    return merged_df


def validate_pacbio_splices(deduplicated_strand_resolved_consensus, db_path, tx_dict, metrics, require_pacbio_support):
    """
    Tag transcripts as having PacBio support.
    If users passed the --require-pacbio-support, remove any transcript which does not have support.
    """
    iso_txs = tools.sqlInterface.load_isoseq_txs(db_path)
    tx_ids, _ = list(zip(*deduplicated_strand_resolved_consensus))
    txs = [tx_dict[tx_id] for tx_id in tx_ids]
    clustered = tools.transcripts.cluster_txs(txs + iso_txs)
    divided_clusters = tools.transcripts.divide_clusters(clustered, tx_ids)
    subset_matches = tools.transcripts.calculate_subset_matches(divided_clusters)
    # invert the subset_matches to extract all validated tx_ids
    validated_ids = set()
    for tx_list in subset_matches.values():
        for tx in tx_list:
            validated_ids.add(tx.name)
    # begin resolving
    pb_resolved_consensus = []
    for tx_id, d in deduplicated_strand_resolved_consensus:
        if tx_id in validated_ids:
            d['pacbio_isoform_supported'] = True
            metrics['IsoSeq Transcript Validation'][True] += 1
            pb_resolved_consensus.append([tx_id, d])
        elif require_pacbio_support is False:
            d['pacbio_isoform_supported'] = False
            metrics['IsoSeq Transcript Validation'][False] += 1
            pb_resolved_consensus.append([tx_id, d])
        # if require_pacbio_support is True, then we don't save this transcript
    return pb_resolved_consensus


def incorporate_tx(best_rows, gene_id, metrics, hints_db_has_rnaseq):
    """incorporate a transcript into the consensus set, storing metrics."""
    best_series = best_rows.iloc[0]
    transcript_modes = evaluate_ties(best_rows)
    # construct the tags for this transcript
    d = {'source_transcript': best_series.TranscriptId,
         'source_transcript_name': best_series.TranscriptName,
         'source_gene': gene_id,
         'score': int(10 * round(best_series.AlnGoodness_mRNA, 3)),
         'transcript_modes': transcript_modes,
         'gene_biotype': best_series.GeneBiotype,
         'transcript_biotype': best_series.TranscriptBiotype,
         'alignment_id': str(best_series.AlignmentId),
         'frameshift': str(best_series.get('Frameshift', None)),
         'exon_annotation_support': ','.join(map(str, best_series.ExonAnnotSupport)),
         'intron_annotation_support': ','.join(map(str, best_series.IntronAnnotSupport)),
         'transcript_class': 'ortholog',
         'valid_start': bool(best_series.ValidStart),
         'valid_stop': bool(best_series.ValidStop),
         'adj_start': best_series.AdjStart_mRNA,
         'adj_stop': best_series.AdjStop_mRNA,
         'proper_orf': bool(best_series.ProperOrf)}
    # incorporate any extra tags
    for key, val in tools.misc.parse_gff_attr_line(best_series.ExtraTags).items():
        d[key] = val
    if hints_db_has_rnaseq is True:
        d['exon_rna_support'] = ','.join(map(str, best_series.ExonRnaSupport))
        d['intron_rna_support'] = ','.join(map(str, best_series.IntronRnaSupport))
    if best_series.Paralogy is not None:
        d['paralogy'] = best_series.Paralogy
    if best_series.UnfilteredParalogy is not None:
        d['unfiltered_paralogy'] = best_series.UnfilteredParalogy
    if best_series.GeneAlternateLoci is not None:
        d['gene_alternate_contigs'] = best_series.GeneAlternateLoci
    if best_series.CollapsedGeneIds is not None:
        d['collapsed_gene_ids'] = best_series.CollapsedGeneIds
    if best_series.CollapsedGeneNames is not None:
        d['collapsed_gene_names'] = best_series.CollapsedGeneNames
    if best_series.PossibleSplitGeneLocations is not None:
        d['possible_split_gene_locations'] = best_series.PossibleSplitGeneLocations
    if best_series.GeneName is not None:
        d['source_gene_common_name'] = best_series.GeneName
    # add information to the overall metrics
    if best_series.TranscriptBiotype == 'protein_coding':
        metrics['Transcript Modes'][transcript_modes] += 1
    metrics['Coverage'][best_series.TranscriptBiotype].append(best_series.AlnCoverage_mRNA)
    metrics['Identity'][best_series.TranscriptBiotype].append(best_series.AlnIdentity_mRNA)
    metrics['Splice Support'][best_series.TranscriptBiotype].append(best_series.IntronRnaSupportPercent)
    metrics['Exon Support'][best_series.TranscriptBiotype].append(best_series.ExonRnaSupportPercent)
    metrics['Splice Annotation Support'][best_series.TranscriptBiotype].append(best_series.IntronAnnotSupportPercent)
    metrics['Exon Annotation Support'][best_series.TranscriptBiotype].append(best_series.ExonAnnotSupportPercent)
    metrics['Original Introns'][best_series.TranscriptBiotype].append(best_series.OriginalIntronsPercent_mRNA)
    return best_series.AlignmentId, d


def evaluate_ties(best_rows):
    """Find out how many transcript modes agreed on this"""
    return ','.join(sorted(set([tools.nameConversions.alignment_type(x) for x in best_rows.AlignmentId])))


def find_novel(db_path, tx_dict, consensus_dict, ref_df, metrics, gene_biotype_map, denovo_num_introns,
               in_species_rna_support_only, denovo_tx_modes, denovo_splice_support, denovo_exon_support,
               denovo_ignore_novel_genes, denovo_novel_end_distance, denovo_allow_unsupported,
               denovo_allow_bad_annot_or_tm, denovo_only_novel_genes, denovo_allow_novel_ends):
    """
    Finds novel loci, builds their attributes. Only calls novel loci if they have sufficient intron and splice support
    as defined by the user.

    Putative novel loci can fall into three categories:
    1) PossibleParlog -- the transcript has no assigned genes, but has alternative genes
    2) PoorMapping -- the transcript has no assigned or alternative genes but has exons/introns supported by annotation
    3) PutativeNovel -- the transcript has no matches to the reference
    4) PossibleFusion -- the transcript was flagged in parent finding as a fusion, and has all valid splices

    If the denovo_ignore_novel_genes flag is set, only incorporate PossibleParalog into the final consensus.

    Also finds novel splice junctions in CGP/PB transcripts. A novel splice junction is defined as a splice which
    homGeneMapping did not map over and which is supported by RNA-seq.
    """
    def is_novel(s):
        """
        Determine if this transcript is possibly novel. If it is assigned a gene ID, pass this off to
         is_novel_supported()
        """
        if s.AssignedGeneId is not None:
            return is_novel_supported(s)
        if denovo_allow_bad_annot_or_tm is False and s.ResolutionMethod == 'badAnnotOrTm':
            return None
        elif s.ResolutionMethod == 'ambiguousOrFusion' and s.IntronRnaSupportPercent != 100:
            return None
        # validate the support level
        intron = s.IntronRnaSupportPercent if in_species_rna_support_only else s.AllSpeciesIntronRnaSupportPercent
        exon = s.ExonRnaSupportPercent if in_species_rna_support_only else s.AllSpeciesExonRnaSupportPercent
        if intron < denovo_splice_support or exon < denovo_exon_support:
            return None
        # if we previously flagged this as ambiguousOrFusion, propagate this tag
        if s.ResolutionMethod == 'ambiguousOrFusion':
            return 'possible_fusion'
        elif s.ResolutionMethod == 'badAnnotOrTm':
            return 'bad_annot_or_tm'
        # if we have alternatives, this is not novel but could be a gene family expansion
        elif s.AlternativeGeneIds is not None:
            return 'possible_paralog'
        # this may be a poor mapping
        elif bool(s.ExonAnnotSupportPercent > 0 or s.CdsAnnotSupportPercent > 0 or s.IntronAnnotSupportPercent > 0):
            return 'poor_alignment'
        # this is looking pretty novel, could still be a mapping problem in a complex region though
        else:
            return 'putative_novel'

    def is_novel_supported(s):
        """Is this CGP/PB transcript with an assigned gene ID supported and have a novel splice?"""
        denovo_tx_obj = tx_dict[s.AlignmentId]
        if len(denovo_tx_obj.intron_intervals) < denovo_num_introns:
            return None
        elif in_species_rna_support_only and s.ExonRnaSupportPercent <= denovo_exon_support or \
                        s.IntronRnaSupportPercent <= denovo_splice_support:
            return None
        elif in_species_rna_support_only is False and s.AllSpeciesExonRnaSupportPercent <= denovo_exon_support or \
                        s.AllSpeciesIntronRnaSupportPercent <= denovo_splice_support:
            return None
        # look for splices that are not supported by the reference annotation
        # these splices may or may not be supported by RNA-seq based on the denovo_allow_unsupported flag
        new_supported_splices = set()
        intron_vector = s.IntronRnaSupport if in_species_rna_support_only else s.AllSpeciesIntronRnaSupport
        for intron, rna in zip(*[denovo_tx_obj.intron_intervals, intron_vector]):
            if intron not in existing_splices:
                if denovo_allow_unsupported is True or rna > 0:
                    new_supported_splices.add(intron)
        if len(new_supported_splices) == 0:
            return None
        # if any splices are both not supported by annotation and supported by RNA, call this as novel
        if any(annot == 0 and i in new_supported_splices for i, annot in zip(*[denovo_tx_obj.intron_intervals,
                                                                               s.IntronAnnotSupport])):
            metrics['Transcript Modes'][tools.nameConversions.alignment_type(s.AlignmentId)] += 1
            tx_class = 'putative_novel_isoform'
        # if any splices are new, and supported by RNA-seq call this poor alignment
        else:
            tx_class = 'poor_alignment'
        return tx_class

    def has_novel_ends(s):
        """Does this transcript have any novel ends we want to retain? Does not apply to augCGP"""
        if tools.nameConversions.aln_id_is_cgp(s.AlignmentId):
            return pd.Series([None, None, s.TranscriptClass])
        denovo_tx_obj = tx_dict[s.AlignmentId]
        five_p = denovo_tx_obj.get_5p_interval()
        three_p = denovo_tx_obj.get_3p_interval()
        five_p_matches = tools.intervals.interval_not_within_wiggle_room_intervals(existing_5p[denovo_tx_obj.chromosome],
                                                                                   five_p, denovo_novel_end_distance)
        three_p_matches = tools.intervals.interval_not_within_wiggle_room_intervals(existing_5p[denovo_tx_obj.chromosome],
                                                                                    three_p, denovo_novel_end_distance)
        if denovo_allow_novel_ends is False:
            tx_class = s.TranscriptClass
        else:
            tx_class = 'putative_novel_isoform' if s.TranscriptClass is None and (five_p_matches or three_p_matches) else s.TranscriptClass
        return pd.Series([five_p_matches, three_p_matches, tx_class])

    denovo_hgm_df = pd.concat([load_hgm_vectors(db_path, tx_mode) for tx_mode in denovo_tx_modes])
    # remove the TranscriptId and GeneId columns so they can be populated by others
    denovo_hgm_df = denovo_hgm_df.drop(['GeneId', 'TranscriptId'], axis=1)
    # load the alignment metrics data
    denovo_alt_names = load_alt_names(db_path, denovo_tx_modes)
    denovo_df = pd.merge(denovo_hgm_df, denovo_alt_names, on='AlignmentId')
    common_name_map = dict(list(zip(ref_df.GeneId, ref_df.GeneName)))

    denovo_df['CommonName'] = [common_name_map.get(x, None) for x in denovo_df.AssignedGeneId]
    denovo_df['GeneBiotype'] = [gene_biotype_map.get(x, None) for x in denovo_df.AssignedGeneId]

    # if we have an external reference, try to incorporate those names as well
    if 'exRef' in denovo_tx_modes:
        def add_exref_ids(s):
            if s.AlignmentId in exref_common_name_map:
                # if we have an assigned gene ID, defer the gene biotype to that but retain transcript biotype
                if s.AssignedGeneId is None:
                    return pd.Series([exref_common_name_map[s.AlignmentId], exref_gene_biotype_map[s.AlignmentId]])
                else:
                    return pd.Series([exref_common_name_map[s.AlignmentId], s.GeneBiotype])
            # pass along the original data
            return pd.Series([s.CommonName, s.GeneBiotype])
        exref_annot = tools.sqlInterface.load_annotation(db_path)
        exref_common_name_map = dict(list(zip(exref_annot.TranscriptId, exref_annot.GeneName)))
        exref_gene_biotype_map = dict(list(zip(exref_annot.TranscriptId, exref_annot.GeneBiotype)))
        denovo_df[['CommonName', 'GeneBiotype']] = denovo_df.apply(add_exref_ids, axis=1)
        exref_annot = exref_annot.set_index('TranscriptId')

    # extract all splices, 5' and 3' ends we have already seen
    existing_splices = set()
    existing_5p = collections.defaultdict(set)
    existing_3p = collections.defaultdict(set)
    for consensus_tx in consensus_dict:
        tx_obj = tx_dict[consensus_tx]
        existing_splices.update(tx_obj.intron_intervals)
        existing_5p[tx_obj.chromosome].add(tx_obj.get_5p_interval())
        existing_3p[tx_obj.chromosome].add(tx_obj.get_3p_interval())

    # apply the novel finding functions
    denovo_df['TranscriptClass'] = denovo_df.apply(is_novel, axis=1)
    denovo_df[['Novel5pCap', 'NovelPolyA', 'TranscriptClass']] = denovo_df.apply(has_novel_ends, axis=1)
    # types of transcripts for later
    denovo_df['TranscriptMode'] = [tools.nameConversions.alignment_type(aln_id) for aln_id in denovo_df.AlignmentId]
    # filter out non-novel as well as fusions
    filtered_denovo_df = denovo_df[(~denovo_df.TranscriptClass.isnull())]
    filtered_denovo_df = filtered_denovo_df[filtered_denovo_df.TranscriptClass != 'possible_fusion']
    # fill in missing fields for novel loci
    filtered_denovo_df['GeneBiotype'] = filtered_denovo_df['GeneBiotype'].fillna('unknown_likely_coding')
    # filter out novel if requested by user
    if denovo_ignore_novel_genes is True:
        filtered_denovo_df = filtered_denovo_df[(filtered_denovo_df.TranscriptClass == 'possible_paralog') |
                                                (filtered_denovo_df.TranscriptClass == 'putative_novel_isoform')]
    elif denovo_only_novel_genes is True:
        filtered_denovo_df = filtered_denovo_df[~((filtered_denovo_df.TranscriptClass == 'possible_paralog') |
                                                (filtered_denovo_df.TranscriptClass == 'putative_novel_isoform'))]

    # construct aln_id -> features map to return
    denovo_tx_dict = {}
    for _, s in filtered_denovo_df.iterrows():
        aln_id = s.AlignmentId
        tx_mode = s.TranscriptMode
        denovo_tx_dict[aln_id] = {'source_gene': s.AssignedGeneId,
                                  'transcript_class': s.TranscriptClass,
                                  'novel_5p_cap': s.Novel5pCap,
                                  'novel_poly_a': s.NovelPolyA,
                                  'transcript_biotype': s.GeneBiotype,
                                  'gene_biotype': s.GeneBiotype,
                                  'intron_rna_support': ','.join(map(str, s.IntronRnaSupport)),
                                  'exon_rna_support': ','.join(map(str, s.ExonRnaSupport)),
                                  'transcript_modes': tx_mode,
                                  'exon_annotation_support': ','.join(map(str, s.ExonAnnotSupport)),
                                  'intron_annotation_support': ','.join(map(str, s.IntronAnnotSupport)),
                                  'alignment_id': aln_id,
                                  'source_gene_common_name': s.CommonName}

        # bring in extra tags for exRef
        if tools.nameConversions.aln_id_is_exref(aln_id):
            tags = exref_annot.loc[aln_id].ExtraTags
            if len(tags) > 0:
                for key, val in tools.misc.parse_gff_attr_line(tags).items():
                    denovo_tx_dict[aln_id][key] = val

        # record some metrics
        metrics['denovo'][tx_mode][s.TranscriptClass.replace('_', ' ').capitalize()] += 1
        metrics['Transcript Modes'][tx_mode] += 1
        metrics['Splice Support']['unknown_likely_coding'].append(s.IntronRnaSupportPercent)
        metrics['Exon Support']['unknown_likely_coding'].append(s.ExonRnaSupportPercent)

    # record how many of each type we threw out
    for tx_mode, df in denovo_df.groupby('TranscriptMode'):
        metrics['denovo'][tx_mode]['Discarded'] = len(df[(~df.TranscriptClass.isnull()) | (~df.Novel5pCap.isnull())
                                                         | (~df.NovelPolyA.isnull())])
    return denovo_tx_dict


def deduplicate_consensus(consensus_dict, tx_dict, metrics):
    """
    In the process of consensus building, we may find that we have ended up with more than one transcript for a gene
    that are actually identical. Remove these, picking the best based on their score, favoring the transcript
    whose biotype matches the parent.
    """
    def resolve_duplicate(tx_list, consensus_dict):
        biotype_txs = [tx for tx in tx_list if
                       consensus_dict[tx].get('gene_biotype', None) == consensus_dict[tx].get('transcript_biotype',
                                                                                              None)]
        if len(biotype_txs) > 0:
            tx_list = biotype_txs
        sorted_scores = sorted([[tx, consensus_dict[tx].get('score', 0)] for tx in tx_list],
                               key=lambda tx_s1: -tx_s1[1])
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
    for tx_list in duplicates.values():
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
    deduplicated_strand_resolved_consensus = []
    for tx_id, attrs in deduplicated_consensus.items():
        tx_obj = tx_dict[tx_id]
        # don't try to resolve novel genes
        source_gene = attrs['source_gene']
        if source_gene is not None:
            gene_dict[source_gene].append([tx_obj, attrs])
        else:
            deduplicated_strand_resolved_consensus.append([tx_obj.name, attrs])

    for gene in gene_dict:
        tx_objs, attrs = list(zip(*gene_dict[gene]))
        if len(set(tx_obj.strand for tx_obj in tx_objs)) > 1:
            strand_scores = collections.Counter()
            for tx_obj, attrs in gene_dict[gene]:
                strand_scores[tx_obj.strand] += attrs.get('score', 0)
            best_strand = sorted(strand_scores.items())[1][0]
            for tx_obj, attrs in gene_dict[gene]:
                if tx_obj.strand == best_strand:
                    deduplicated_strand_resolved_consensus.append([tx_obj.name, attrs])
                else:
                    metrics['Discarded by strand resolution'] += 1
        else:
            deduplicated_strand_resolved_consensus.extend([[tx_obj.name, attrs] for tx_obj, attrs in gene_dict[gene]])
    return deduplicated_strand_resolved_consensus


###
# This section involves final filtering for overlapping genes caused by incorporating predictions
###


def resolve_overlapping_cds_intervals(overlapping_gene_distance, deduplicated_strand_resolved_consensus, tx_dict):
    """
    Resolves overlapping genes that are the result of integrating gene predictions
    """
    # first, write genePred
    attr_df = []
    with tools.fileOps.TemporaryFilePath() as tmp_gp, tools.fileOps.TemporaryFilePath() as tmp_clustered:
        with open(tmp_gp, 'w') as outf:
            for tx_id, attrs in deduplicated_strand_resolved_consensus:
                tx_obj = tx_dict[tx_id]
                tools.fileOps.print_row(outf, tx_obj.get_gene_pred())
                attr_df.append([tx_id, attrs['transcript_class'], attrs['gene_biotype'],
                                attrs.get('source_gene', tx_obj.name2), attrs.get('score', None)])
        # cluster
        cmd = ['clusterGenes', '-cds', f'-minOverlappingBases={overlapping_gene_distance}',
               tmp_clustered, 'no', tmp_gp]
        tools.procOps.run_proc(cmd)
        cluster_df = pd.read_csv(tmp_clustered, sep='\t')
    attr_df = pd.DataFrame(attr_df, columns=['transcript_id', 'transcript_class', 'gene_biotype', 'gene_id', 'score'])
    m = attr_df.merge(cluster_df, left_on='transcript_id', right_on='gene')  # gene is transcript ID

    to_remove = set()  # list of transcript IDs to remove
    for cluster_id, group in m.groupby('#cluster'):
        if len(set(group['gene_id'])) > 1:
            if 'unknown_likely_coding' in set(group['gene_biotype']):  # pick longest ORF
                orfs = {tx_id: tx_dict[tx_id].cds_size for tx_id in group['transcript_id']}
                best_tx = sorted(iter(orfs.items()), key=lambda x: x[1])[-1][0]
                tx_df = group[group.transcript_id == best_tx].iloc[0]
                best_gene = tx_df.gene_id
            else:  # pick highest average score
                avg_scores = group[['gene_id', 'score']].groupby('gene_id', as_index=False).mean()
                best_gene = avg_scores.sort_values('score', ascending=False).iloc[0]['gene_id']
            to_remove.update(set(group[group.gene_id != best_gene].transcript_id))

    return [[tx_id, attrs] for tx_id, attrs in deduplicated_strand_resolved_consensus if tx_id not in to_remove]


###
# Metrics calculations on final set
###


def calculate_completeness(final_consensus, metrics):
    """calculates final completeness to make arithmetic easier"""
    genes = collections.defaultdict(set)
    txs = collections.Counter()
    for aln_id, c in final_consensus:
        # don't count novel transcripts towards completeness
        if tools.nameConversions.aln_id_is_cgp(aln_id) or tools.nameConversions.aln_id_is_pb(aln_id):
            continue
        genes[c['gene_biotype']].add(c['source_gene'])
        txs[c['transcript_biotype']] += 1
    genes = {biotype: len(gene_list) for biotype, gene_list in genes.items()}
    metrics['Completeness'] = {'Gene': genes, 'Transcript': txs}


def calculate_improvement_metrics(final_consensus, scored_df, tm_eval_df, hgm_df, metrics):
    """For coding transcripts, how much did we improve the metrics?"""
    tm_df = tm_eval_df.reset_index()[['TransMapOriginalIntronsPercent', 'TranscriptId']]
    hgm_df_subset = hgm_df[hgm_df['AlignmentId'].apply(tools.nameConversions.aln_id_is_transmap)]
    hgm_df_subset = hgm_df_subset[['TranscriptId', 'IntronAnnotSupportPercent', 'IntronRnaSupportPercent']]
    tm_df = pd.merge(tm_df, hgm_df_subset, on='TranscriptId')
    df = pd.merge(tm_df, scored_df.reset_index(), on='TranscriptId', suffixes=['TransMap', ''])
    df = df.drop_duplicates(subset='AlignmentId')  # why do I need to do this?
    df = df.set_index('AlignmentId')
    metrics['Evaluation Improvement'] = {'changes': [], 'unchanged': 0}
    for aln_id, c in final_consensus:
        if c['transcript_biotype'] != 'protein_coding':
            continue
        elif 'exRef' in c['transcript_modes'] or 'augPB' in c['transcript_modes'] or 'augCGP' in c['transcript_modes']:
            continue
        elif 'transMap' in c['transcript_modes']:
            metrics['Evaluation Improvement']['unchanged'] += 1
            continue
        tx_s = df.loc[aln_id]
        metrics['Evaluation Improvement']['changes'].append([tx_s.TransMapOriginalIntronsPercent,
                                                             tx_s.IntronAnnotSupportPercentTransMap,
                                                             tx_s.IntronRnaSupportPercentTransMap,
                                                             tx_s.OriginalIntronsPercent_mRNA,
                                                             tx_s.IntronAnnotSupportPercent,
                                                             tx_s.IntronRnaSupportPercent,
                                                             tx_s.TransMapGoodness,
                                                             tx_s.AlnGoodness_mRNA])


def calculate_indel_metrics(final_consensus, eval_df, metrics):
    """How many transcripts in the final consensus have indels? How many did we have in transMap?"""
    if len(eval_df) == 0:  # edge case where no transcripts hit indel filters
        metrics['transMap Indels'] = {}
        metrics ['Consensus Indels'] = {}
        return
    eval_df_transmap = eval_df[eval_df['AlignmentId'].apply(tools.nameConversions.aln_id_is_transmap)]
    tm_vals = eval_df_transmap.set_index('AlignmentId').sum(axis=0)
    tm_vals = 100.0 * tm_vals / len(set(eval_df_transmap.index))
    metrics['transMap Indels'] = tm_vals.to_dict()
    consensus_ids = set(list(zip(*final_consensus))[0])
    consensus_vals = eval_df[eval_df['AlignmentId'].isin(consensus_ids)].set_index('AlignmentId').sum(axis=0)
    consensus_vals = 100.0 * consensus_vals / len(final_consensus)
    metrics['Consensus Indels'] = consensus_vals.to_dict()


###
# Writing output genePred and GFF3
###


def write_consensus_gps(consensus_gp, consensus_gp_info, final_consensus, tx_dict, genome, gene_offset=0):
    """
    Write the resulting gp + gp_info, generating genome-specific unique identifiers
    """
    # keeps track of gene # -> ID mappings
    genes_seen = collections.defaultdict(dict)
    gene_count = gene_offset
    consensus_gene_dict = DefaultOrderedDict(lambda: DefaultOrderedDict(list))  # used to make gff3 next
    gp_infos = []
    consensus_gp_target = luigi.LocalTarget(consensus_gp)
    with consensus_gp_target.open('w') as out_gp:
        for tx_count, (tx, attrs) in enumerate(final_consensus, 1):
            attrs = attrs.copy()
            tx_obj = tx_dict[tx]
            name = id_template.format(genome=genome, tag_type='T', unique_id=tx_count)
            score = int(round(attrs.get('score', 0)))
            source_gene = attrs['source_gene']
            if source_gene is None:
                source_gene = tx_obj.name2
            if source_gene not in genes_seen[tx_obj.chromosome]:
                gene_count += 1
                genes_seen[tx_obj.chromosome][source_gene] = gene_count
            gene_id = genes_seen[tx_obj.chromosome][source_gene]
            name2 = id_template.format(genome=genome, tag_type='G', unique_id=gene_id)
            out_gp.write('\t'.join(tx_obj.get_gene_pred(name=name, name2=name2, score=score)) + '\n')
            attrs['transcript_id'] = name
            attrs['gene_id'] = name2
            gp_infos.append(attrs)
            consensus_gene_dict[tx_obj.chromosome][name2].append([tx_obj, attrs])
    gp_info_df = pd.DataFrame(gp_infos)
    gp_info_df = gp_info_df.set_index(['gene_id', 'transcript_id'])
    # its possible alternative_source_transcripts did not end up in the final result, so add it
    if 'alternative_source_transcripts' not in gp_info_df.columns:
        gp_info_df['alternative_source_transcripts'] = ['N/A'] * len(gp_info_df)
    with luigi.LocalTarget(consensus_gp_info).open('w') as outf:
        gp_info_df.to_csv(outf, sep='\t', na_rep='N/A')
    return consensus_gene_dict


def write_consensus_gff3(consensus_gene_dict, consensus_gff3):
    """
    Write the consensus set in gff3 format
    """
    def convert_attrs(attrs, id_field):
        """converts the attrs dict to a attributes field. assigns name to the gene common name for display"""
        attrs['ID'] = id_field
        if 'score' in attrs:
            score = 10 * attrs['score']
            del attrs['score']
        else:
            score = '.'
        if 'source_gene_common_name' in attrs and isinstance(attrs['source_gene_common_name'], str):
            attrs['Name'] = attrs['source_gene_common_name']
        else:
            attrs['Name'] = attrs['gene_id']
        # convert empty strings into nan
        attrs_str = []
        for key, val in attrs.items():
            val = str(val)
            if len(val) == 0:
                val = 'nan'
            val = str(val).replace('=', '%3D').replace(';', '%3B')
            key = key.replace('=', '%3D').replace(';', '%3B')
            attrs_str.append(f"{key}={val}")
        return score, ";".join(attrs_str)

    def find_feature_support(attrs, feature, i):
        """Extracts the boolean value from the comma delimited string"""
        try:
            vals = list(map(bool, attrs[feature].split(',')))
        except KeyError:
            return 'N/A'
        return vals[i]

    def generate_gene_record(chrom, tx_objs, gene_id, attrs_list):
        """calculates the gene interval for this list of tx"""
        def find_all_tx_modes(attrs_list):
            tx_modes = set()
            for attrs in attrs_list:
                tx_modes.update(attrs['transcript_modes'].split(','))
            return ','.join(tx_modes)

        intervals = set()
        for tx in tx_objs:
            intervals.update(tx.exon_intervals)
        intervals = sorted(intervals)
        strand = tx_objs[0].strand
        # subset the attrs to gene fields
        attrs = attrs_list[0]
        useful_keys = ['source_gene_common_name', 'source_gene', 'gene_biotype',
                       'alternative_source_transcripts', 'gene_alternate_contigs',
                       'gene_name', 'gene_id']
        attrs = {key: attrs[key] for key in useful_keys if key in attrs}
        attrs['transcript_modes'] = find_all_tx_modes(attrs_list)
        score, attrs_field = convert_attrs(attrs, gene_id)
        return [chrom, 'CAT', 'gene', intervals[0].start + 1, intervals[-1].stop, score, strand, '.', attrs_field]

    def generate_transcript_record(chrom, tx_obj, attrs):
        """generates transcript records, calls generate_exon_records to generate those too"""
        tx_id = attrs['transcript_id']
        attrs['Parent'] = attrs['gene_id']
        score, attrs_field = convert_attrs(attrs, tx_id)
        yield [chrom, 'CAT', 'transcript', tx_obj.start + 1, tx_obj.stop, score, tx_obj.strand, '.', attrs_field]
        # hack to remove the frameshift field from lower objects
        # TODO: record the actual exon with the frameshift.
        if 'frameshift' in attrs:
            del attrs['frameshift']
        for line in generate_intron_exon_records(chrom, tx_obj, tx_id, attrs):
            yield line
        if tx_obj.cds_size > 3:
            for line in generate_start_stop_codon_records(chrom, tx_obj, tx_id, attrs):
                yield line

    def generate_intron_exon_records(chrom, tx_obj, tx_id, attrs):
        """generates intron and exon records"""
        attrs['Parent'] = tx_id
        # exon records
        cds_i = 0  # keep track of position of CDS in case of entirely non-coding exons
        for i, (exon, exon_frame) in enumerate(zip(*[tx_obj.exon_intervals, tx_obj.exon_frames])):
            attrs['rna_support'] = find_feature_support(attrs, 'exon_rna_support', i)
            attrs['reference_support'] = find_feature_support(attrs, 'exon_annotation_support', i)
            score, attrs_field = convert_attrs(attrs, 'exon:{}:{}'.format(tx_id, i))
            yield [chrom, 'CAT', 'exon', exon.start + 1, exon.stop, score, exon.strand, '.', attrs_field]
            cds_interval = exon.intersection(tx_obj.coding_interval)
            if cds_interval is not None:
                score, attrs_field = convert_attrs(attrs, 'CDS:{}:{}'.format(tx_id, cds_i))
                cds_i += 1
                yield [chrom, 'CAT', 'CDS', cds_interval.start + 1, cds_interval.stop, score, exon.strand,
                       tools.transcripts.convert_frame(exon_frame), attrs_field]

        # intron records
        for i, intron in enumerate(tx_obj.intron_intervals):
            if len(intron) == 0:
                continue
            attrs['rna_support'] = find_feature_support(attrs, 'intron_rna_support', i)
            attrs['reference_support'] = find_feature_support(attrs, 'intron_annotation_support', i)
            score, attrs_field = convert_attrs(attrs, 'intron:{}:{}'.format(tx_id, i))
            yield [chrom, 'CAT', 'intron', intron.start + 1, intron.stop, score, intron.strand, '.', attrs_field]

    def generate_start_stop_codon_records(chrom, tx_obj, tx_id, attrs):
        """generate start/stop codon GFF3 records, handling frame appropriately"""
        if attrs.get('valid_start') is True:
            score, attrs_field = convert_attrs(attrs, 'start_codon:{}'.format(tx_id))
            for interval in tx_obj.get_start_intervals():
                yield [chrom, 'CAT', 'start_codon', interval.start + 1, interval.stop, score, tx_obj.strand,
                       interval.data, attrs_field]
        if attrs.get('valid_stop') is True:
            score, attrs_field = convert_attrs(attrs, 'stop_codon:{}'.format(tx_id))
            for interval in tx_obj.get_stop_intervals():
                yield [chrom, 'CAT', 'stop_codon', interval.start + 1, interval.stop, score, tx_obj.strand,
                       interval.data, attrs_field]

    # main gff3 writing logic
    consensus_gff3 = luigi.LocalTarget(consensus_gff3)
    with consensus_gff3.open('w') as out_gff3:
        out_gff3.write('##gff-version 3\n')
        for chrom in sorted(consensus_gene_dict):
            for gene_id, tx_list in consensus_gene_dict[chrom].items():
                tx_objs, attrs_list = list(zip(*tx_list))
                tx_lines = [generate_gene_record(chrom, tx_objs, gene_id, attrs_list)]
                for tx_obj, attrs in tx_list:
                    tx_lines.extend(list(generate_transcript_record(chrom, tx_obj, attrs)))
                tx_lines = sorted(tx_lines, key=lambda l: l[3])
                tools.fileOps.print_rows(out_gff3, tx_lines)


def write_consensus_fastas(consensus_gene_dict, consensus_fasta, consensus_protein_fasta, fasta):
    """Write FASTA records for both"""
    seq_dict = tools.bio.get_sequence_dict(fasta)
    consensus_fasta = luigi.LocalTarget(consensus_fasta)
    consensus_protein_fasta = luigi.LocalTarget(consensus_protein_fasta)
    with consensus_fasta.open('w') as cfa, consensus_protein_fasta.open('w') as cpfa:
        for chrom in sorted(consensus_gene_dict):
            for gene_id, tx_list in consensus_gene_dict[chrom].items():
                for tx_obj, attrs in tx_list:
                    tools.bio.write_fasta(cfa, attrs['transcript_id'], tx_obj.get_mrna(seq_dict))
                    if tx_obj.cds_size > 0:
                        tools.bio.write_fasta(cpfa, attrs['transcript_id'], tx_obj.get_protein_sequence(seq_dict))
