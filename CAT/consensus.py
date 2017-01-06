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
5. alternative_source_transcripts: A comma separated list of alternate IDs for this transcript
6. failed_gene: This transcript is the single representative for a failed transcript
7. transcript_class: One of possible_paralog, poor_alignment, putative_novel, putative_novel_isoform, failing, passing
8. paralogy: The number of paralogs that were mapped over when transMap mapped this transcript
9. gene_biotype: gene biotype
10. transcript_biotype: transcript biotype
11. paralog_status: confident if the paralog was confidently resolved, not_confident if it was not
12. alternative_source_transcripts: Other possible transcripts, if this was collapsed as the result of deduplication
13. gene_alternate_contigs: If the --resolve-split-genes flag was set, contigs that this gene was also found on are
    comma separated in this tag.
14: transcript_modes: The mode(s) that generated this transcript
"""
import collections
import copy
import luigi
import logging
import pandas as pd

import tools.intervals
import tools.mathOps
import tools.fileOps
import tools.sqlInterface
import tools.transcripts
import tools.nameConversions
from tools.defaultOrderedDict import DefaultOrderedDict

logger = logging.getLogger(__name__)

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
    # gene transcript map to iterate over so that we capture missing gene information
    gene_transcript_map = tools.sqlInterface.get_gene_transcript_map(args.ref_db_path)
    gene_biotype_map = tools.sqlInterface.get_gene_biotype_map(args.ref_db_path)
    transcript_biotype_map = tools.sqlInterface.get_transcript_biotype_map(args.ref_db_path)
    common_name_map = dict(zip(*[ref_df.GeneId, ref_df.GeneName]))
    # load transMap evaluation data
    tm_eval_df = load_transmap_evals(args.db_path)
    # load the transMap percent ID cutoff previously identified
    coding_cutoff = tools.sqlInterface.load_tm_fit(args.db_path)
    # load the homGeneMapping data for transMap/augTM/augTMR
    tx_modes = [x for x in args.tx_modes if x in ['transMap', 'augTM', 'augTMR']]
    hgm_df = pd.concat([load_hgm_vectors(args.db_path, tx_mode) for tx_mode in tx_modes])
    # load the alignment metrics data
    metrics_df = pd.concat([load_metrics_from_db(args.db_path, tx_mode, 'mRNA') for tx_mode in tx_modes])
    # combine and filter the results
    coding_df, non_coding_df = combine_and_filter_dfs(hgm_df, metrics_df, tm_eval_df, ref_df,
                                                      args.intron_rnaseq_support, args.exon_rnaseq_support,
                                                      args.intron_annot_support, args.exon_annot_support,
                                                      args.original_intron_support, coding_cutoff,
                                                      args.minimum_coverage,
                                                      args.in_species_rna_support_only)
    if len(coding_df) + len(non_coding_df) == 0:
        raise RuntimeError('No transcripts pass filtering for species {}. '
                           'Consider lowering requirements. Please see the manual.'.format(args.genome))
    elif len(coding_df) == 0:
        logger.warning('No protein coding transcripts pass filtering for species {}. '
                       'Consider lowering requirements. Please see the manual.'.format(args.genome))
    elif len(non_coding_df) == 0:
        logger.warning('No non-coding transcripts pass filtering for species {}. '
                       'Consider lowering requirements. Please see the manual.'.format(args.genome))
    scored_coding_df, scored_non_coding_df = score_filtered_dfs(coding_df, non_coding_df,
                                                                args.in_species_rna_support_only)
    scored_df = merge_scored_dfs(scored_coding_df, scored_non_coding_df)

    # store some metrics for plotting
    metrics = {'Gene Failed': collections.Counter(),
               'Transcript Failed': collections.Counter(),
               'Transcript Missing': collections.Counter(),
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
               'IsoSeq Transcript Valdiation': collections.Counter()}

    # stores a mapping of alignment IDs to tags for the final consensus set
    consensus_dict = {}

    # if we ran in either denovo mode, load those data and detect novel genes
    if len(args.denovo_tx_modes) > 0:
        metrics['denovo'] = {}
        for tx_mode in args.denovo_tx_modes:
            metrics['denovo'][tx_mode] = {'Possible paralog': 0, 'Poor mapping': 0, 'Putative novel': 0,
                                          'Discarded': 0, 'Novel isoforms': 0}  # isoforms populated later
        denovo_hgm_df = pd.concat([load_hgm_vectors(args.db_path, tx_mode) for tx_mode in args.denovo_tx_modes])
        # remove the TranscriptId and GeneId columns so they can be populated by others
        denovo_hgm_df = denovo_hgm_df.drop(['GeneId', 'TranscriptId'], axis=1)
        # load the alignment metrics data
        denovo_alt_names = load_alt_names(args.db_path, args.denovo_tx_modes)
        denovo_df = pd.merge(denovo_hgm_df, denovo_alt_names, on='AlignmentId')
        find_novel_transcripts(denovo_df, tx_dict, args.denovo_num_introns, args.denovo_splice_support,
                               args.denovo_exon_support, metrics, consensus_dict)
        denovo_df = denovo_df.set_index('AssignedGeneId')

    # main consensus finding logic. Start iterating over each gene then each transcript
    for gene_id, tx_list in gene_transcript_map.iteritems():
        gene_consensus_dict = {}
        gene_df = slice_df(scored_df, gene_id)
        gene_biotype = gene_biotype_map[gene_id]
        # if we have no transcripts, record this gene and all tx's as missing and continue
        if len(gene_df) == 0:
            metrics['Gene Missing'][gene_biotype] += 1
            for tx_id in tx_list:
                tx_biotype = transcript_biotype_map[tx_id]
                metrics['Transcript Missing'][tx_biotype] += 1
            continue
        # evaluate if this gene is failing.
        failed_gene = is_failed_df(gene_df)
        if failed_gene is True:
            aln_id, d = rescue_failed_gene(gene_df, gene_id, metrics, args.hints_db_has_rnaseq)
            gene_consensus_dict[aln_id] = d
            metrics['Gene Failed'][gene_biotype] += 1
        else:  # begin consensus finding for each transcript
            for tx_id in tx_list:
                tx_biotype = transcript_biotype_map[tx_id]
                tx_df = slice_df(gene_df, tx_id)
                if len(tx_df) == 0:  # keep track of isoforms that did not map over
                    metrics['Transcript Missing'][tx_biotype] += 1
                elif is_failed_df(tx_df):  # failed transcripts do not get incorporated
                    metrics['Transcript Failed'][tx_biotype] += 1
                else:
                    best_rows = find_best_score(tx_df)
                    aln_id, d = incorporate_tx(best_rows, gene_id, metrics, args.hints_db_has_rnaseq, failed_gene=False)
                    gene_consensus_dict[aln_id] = d
        if len(args.denovo_tx_modes) > 0:
            denovo_gene_df = slice_df(denovo_df, gene_id)
            if len(denovo_gene_df) == 0:
                continue
            denovo_gene_df = denovo_gene_df.set_index('AlignmentId')
            gene_consensus_dict.update(find_novel_splices(gene_consensus_dict, denovo_gene_df, tx_dict, gene_id,
                                                          common_name_map, metrics, failed_gene, gene_biotype,
                                                          args.denovo_num_introns))
        consensus_dict.update(gene_consensus_dict)

    # perform final filtering steps
    deduplicated_consensus = deduplicate_consensus(consensus_dict, tx_dict, metrics)
    deduplicated_strand_resolved_consensus = resolve_opposite_strand(deduplicated_consensus, tx_dict, metrics)

    if 'augPB' in args.denovo_tx_modes:
        deduplicated_strand_resolved_consensus = validate_pacbio_splices(deduplicated_strand_resolved_consensus,
                                                                         args.db_path, tx_dict, metrics,
                                                                         args.require_pacbio_support)

    # sort by genomic interval for prettily increasing numbers
    final_consensus = sorted(deduplicated_strand_resolved_consensus,
                             key=lambda (tx, attrs): (tx_dict[tx].chromosome, tx_dict[tx].start))

    # calculate final gene set completeness
    calculate_completeness(final_consensus, metrics)
    # add some interesting metrics on how much using Augustus modes improved our results
    if 'augTM' or 'augTMR' in tx_modes:
        calculate_improvement_metrics(final_consensus, scored_df, tm_eval_df, metrics)

    # write out results. consensus tx dict has the unique names
    consensus_gene_dict = write_consensus_gps(args.consensus_gp, args.consensus_gp_info,
                                              final_consensus, tx_dict, args.genome)
    write_consensus_gff3(consensus_gene_dict, args.consensus_gff3)

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
    merged = pd.merge(tm_eval, tm_filter_eval, on=['TranscriptId', 'AlignmentId'])
    # remove the AlignmentId column to make the future merges work
    return merged.drop('AlignmentId', axis=1)


def parse_text_vector(s):
    """Used by load_hgm_vectors() and load_metrics() to parse text vectors into numeric lists"""
    if s is None or s == '':
        return []
    elif isinstance(s, float):
        return [s]
    else:
        return map(int, s.split(','))


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
        hgm_df[col] = hgm_df[col].apply(parse_text_vector)
        hgm_df[col + 'Percent'] = hgm_df[col].apply(calculate_vector_support, resolve_nan=1)
    return hgm_df


def load_metrics_from_db(db_path, tx_mode, aln_mode):
    """
    Loads the alignment metrics for the mRNA alignments of transMap/AugustusTM/TMR
    """
    def aggfunc(s):
        """used to aggregate columns. Attempts to convert each cell to a float if possible"""
        try:
            return float(s)
        except ValueError:
            return s.iloc[0]
        except TypeError:
            assert s.iloc[0] is None
            return None

    session = tools.sqlInterface.start_session(db_path)
    metrics_table = tools.sqlInterface.tables[aln_mode][tx_mode]['metrics']
    metrics_df = tools.sqlInterface.load_metrics(metrics_table, session)
    metrics_df = pd.pivot_table(metrics_df, index=['GeneId', 'TranscriptId', 'AlignmentId'], columns='classifier',
                                values='value', fill_value=None, aggfunc=aggfunc).reset_index()
    metrics_df['OriginalIntrons'] = metrics_df['OriginalIntrons'].apply(parse_text_vector)
    metrics_df['OriginalIntronsPercent'] = metrics_df['OriginalIntrons'].apply(calculate_vector_support)
    session.close()
    return metrics_df


def load_alt_names(db_path, denovo_tx_modes):
    """Load the alternative tx tables for augCGP/augPB"""
    session = tools.sqlInterface.start_session(db_path)
    r = []
    for tx_mode in denovo_tx_modes:
        table = tools.sqlInterface.AugCgpAlternativeGenes if tx_mode == 'augCGP' else tools.sqlInterface.AugPbAlternativeGenes
        r.append(tools.sqlInterface.load_alternatives(table, session))
    df = pd.concat(r)
    # rename TranscriptId to AlignmentId. This is all super confusing and silly
    # the reason is that homGeneMapping has the gene -> tx -> aln ID hierarchy we inherit to simplify things
    df.columns = [x if x != 'TranscriptId' else 'AlignmentId' for x in df.columns]
    return df


def combine_and_filter_dfs(hgm_df, metrics_df, tm_eval_df, ref_df, intron_rnaseq_support, exon_rnaseq_support,
                           intron_annot_support, exon_annot_support, original_intron_support, coding_cutoff,
                           minimum_coverage, in_species_rna_support_only):
    """
    Updates the DataFrame based on support levels. Filters based on user-tunable flags for support levels.
    :param hgm_df: df produced by load_hgm_vectors() (all transcripts)
    :param metrics_df: df produced by load_metrics() (coding transcripts only)
    :param tm_eval_df: df produced by load_transmap_evals() (all transcripts)
    :param ref_df: df produced by tools.sqlInterface.load_annotation()
    :param intron_rnaseq_support: Value 0-100. Percent of introns that must be supported by RNAseq
    :param exon_rnaseq_support: Value 0-100. Percent of exons supported by RNA-seq.
    :param intron_annot_support: Value 0-100. Percent of introns supported by the reference.
    :param exon_annot_support: Value 0-100. Percent of exons supported by the reference.
    :param original_intron_support: Value 0-100. Percent of introns that must be supported by this specific annotation.
    :param coding_cutoff: Lognormal fit derived cutoff for coding transcripts. Used to filter TM/TMR alignments.
    :param minimum_coverage: Minimum alignment coverage to be considered.
    :param in_species_rna_support_only: Should we use the homGeneMapping vectors within-species or all-species?
    :return: filtered and merged dataframe
    """
    # remove the start/stop codon information from the ref_df because we don't currently use it and it makes life hard
    ref_df = ref_df.drop(['StartCodon', 'StopCodon'], axis=1)
    # add the reference information to gain biotype information
    hgm_ref_df = pd.merge(hgm_df, ref_df, on=['GeneId', 'TranscriptId'])
    # combine in homGeneMapping results
    hgm_ref_tm_df = pd.merge(hgm_ref_df, tm_eval_df, on=['GeneId', 'TranscriptId'])
    # split merged_df into coding and noncoding
    coding_df = hgm_ref_tm_df[hgm_ref_tm_df.TranscriptBiotype == 'protein_coding']
    non_coding_df = hgm_ref_tm_df[hgm_ref_tm_df.TranscriptBiotype != 'protein_coding']
    # add metrics information to coding df
    coding_df = pd.merge(coding_df, metrics_df, on=['GeneId', 'TranscriptId', 'AlignmentId'])

    # huge ugly filtering expression for coding transcripts
    if in_species_rna_support_only is True:
        filt = ((coding_df.AlnIdentity >= coding_cutoff) &
                (coding_df.AlnCoverage >= minimum_coverage) &
                (coding_df.OriginalIntronsPercent >= original_intron_support) &
                (coding_df.IntronAnnotSupportPercent >= intron_annot_support) &
                (coding_df.IntronRnaSupportPercent >= intron_rnaseq_support) &
                (coding_df.ExonAnnotSupportPercent >= exon_annot_support) &
                (coding_df.ExonRnaSupportPercent >= exon_rnaseq_support))
    else:
        filt = ((coding_df.AlnIdentity >= coding_cutoff) &
                (coding_df.AlnCoverage >= minimum_coverage) &
                (coding_df.OriginalIntronsPercent >= original_intron_support) &
                (coding_df.IntronAnnotSupportPercent >= intron_annot_support) &
                (coding_df.AllSpeciesIntronRnaSupportPercent >= intron_rnaseq_support) &
                (coding_df.ExonAnnotSupportPercent >= exon_annot_support) &
                (coding_df.AllSpeciesExonRnaSupportPercent >= exon_rnaseq_support))
    coding_df = coding_df[filt]

    # huge ugly filtering expression for non coding transcripts
    if in_species_rna_support_only is True:
        filt = ((non_coding_df.TransMapCoverage >= minimum_coverage) &
                (non_coding_df.TransMapOriginalIntronsPercent >= original_intron_support) &
                (non_coding_df.IntronAnnotSupportPercent >= intron_annot_support) &
                (non_coding_df.IntronRnaSupportPercent >= intron_rnaseq_support) &
                (non_coding_df.ExonAnnotSupportPercent >= exon_annot_support) &
                (non_coding_df.ExonRnaSupportPercent >= exon_rnaseq_support))
    else:
        filt = ((non_coding_df.TransMapCoverage >= minimum_coverage) &
                (non_coding_df.TransMapOriginalIntronsPercent >= original_intron_support) &
                (non_coding_df.IntronAnnotSupportPercent >= intron_annot_support) &
                (non_coding_df.AllSpeciesIntronRnaSupportPercent >= intron_rnaseq_support) &
                (non_coding_df.ExonAnnotSupportPercent >= exon_annot_support) &
                (non_coding_df.AllSpeciesExonRnaSupportPercent >= exon_rnaseq_support))
    non_coding_df = non_coding_df[filt]

    return coding_df, non_coding_df


def score_filtered_dfs(coding_df, non_coding_df, in_species_rna_support_only):
    """
    Scores the alignments. The score is the additive combination of the following features:
    1) Alignment goodness.
    2) Intron annotation support.
    3) Exon annotation support.
    4) Original intron support.
    If we have RNA-seq data, the following fields are also incorporated:
    5) Intron RNA-seq support.
    6) Exon RNA-seq support.

    In some future world these features could be combined differently, maybe with some fancy machine learning.

    Returns the dataframe sorted by scores after indexing.
    """
    def score(s):
        goodness = s.AlnGoodness if s.TranscriptBiotype == 'protein_coding' else s.TransMapGoodness
        orig_intron = s.OriginalIntronsPercent if s.TranscriptBiotype == 'protein_coding' else s.TransMapOriginalIntronsPercent
        if in_species_rna_support_only:
            rna_support = s.ExonRnaSupportPercent + s.IntronRnaSupportPercent
        else:
            rna_support = s.AllSpeciesExonRnaSupportPercent + s.AllSpeciesIntronRnaSupportPercent
        return goodness + s.IntronAnnotSupportPercent + s.ExonAnnotSupportPercent + orig_intron + rna_support

    coding_df['TranscriptScore'] = coding_df.apply(score, axis=1)
    non_coding_df['TranscriptScore'] = non_coding_df.apply(score, axis=1)
    return coding_df, non_coding_df


def merge_scored_dfs(scored_coding_df, scored_non_coding_df):
    """Merges the scored dataframes by changing some names around"""
    # for every non-coding TransMap metric, copy it to the other name
    for m in ['Coverage', 'Identity', 'Goodness']:
        scored_non_coding_df['Aln' + m] = scored_non_coding_df['TransMap' + m]
    merged_df = pd.concat([scored_non_coding_df, scored_coding_df])
    merged_df = merged_df.set_index(['GeneId', 'TranscriptId'])
    merged_df = merged_df.sort_values('TranscriptScore', ascending=False)
    return merged_df


def find_novel_transcripts(denovo_df, tx_dict, denovo_num_introns, denovo_splice_support, denovo_exon_support, metrics,
                           consensus_dict):
    """
    Finds novel loci, builds their attributes. Only calls novel loci if they have sufficient intron and splice support
    as defined by the user.

    Putative novel loci can fall into three categories:
    1) PossibleParlog -- the transcript has no assigned genes, but has alternative genes
    2) PoorMapping -- the transcript has no assigned or alternative genes but has exons/introns supported by annotation
    1) PutativeNovel -- the transcript has no matches to the reference

    """
    def is_supported(s, tx):
        """validate support levels for this putative novel transcript"""
        return s.IntronRnaSupportPercent >= denovo_splice_support and \
               s.ExonRnaSupportPercent >= denovo_exon_support and \
               len(tx.intron_intervals) >= denovo_num_introns

    def is_possible_paralog(s):
        """if we have alternative gene IDs this is a possible paralog"""
        return s.AlternativeGeneIds is not None

    def is_poor_alignment(s):
        """If we have no alternative GeneIds, but we have annotation support, this may be a poorly mapped gene"""
        return s.ExonAnnotSupportPercent == 0 and s.CdsAnnotSupportPercent == 0 and s.IntronAnnotSupportPercent == 0

    # novel loci will have None in the AssignedGeneId field but may be non-None in the AlternativeGeneIds field
    for aln_id, df in denovo_df.groupby('AlignmentId'):
        if df.iloc[0].AssignedGeneId is None:
            assert len(df) == 1
            s = df.iloc[0]
            tx = tx_dict[s.AlignmentId]
            tx_mode = tools.nameConversions.alignment_type(aln_id)
            # validate the support level
            if is_supported(s, tx) is False:
                metrics['discarded'][tx_mode]['Discarded'] += 1
            d = {'gene_biotype': 'unknown_likely_coding', 'transcript_biotype': 'unknown_likely_coding'}
            # if we have alternatives, this is not novel but could be a gene family expansion
            if is_possible_paralog(s):
                d['transcript_class'] = 'possible_paralog'
                metrics['denovo'][tx_mode]['Possible paralog'] += 1
                metrics['Transcript Modes'][tx_mode] += 1
            # if we have no alternatives assigned, but we have any sign of mapped over annotations,
            # this may be a poor mapping
            elif is_poor_alignment(s):
                d['transcript_class'] = 'poor_alignment'
                metrics['denovo'][tx_mode]['Poor mapping'] += 1
                metrics['Transcript Modes'][tx_mode] += 1
            # this is looking pretty novel, could still be a mapping problem in a complex region though
            else:
                d['transcript_class'] = 'putative_novel'
                metrics['denovo'][tx_mode]['Putative novel'] += 1
                metrics['Transcript Modes'][tx_mode] += 1
            d['transcript_modes'] = tx_mode
            consensus_dict[aln_id] = d
            metrics['Transcript Modes'][tx_mode] += 1
            d['exon_rna_support'] = ','.join(map(str, s.ExonRnaSupport))
            d['intron_rna_support'] = ','.join(map(str, s.IntronRnaSupport))
            d['exon_annotation_support'] = ','.join(map(str, s.ExonAnnotSupport))
            d['cds_annotation_support'] = ','.join(map(str, s.CdsAnnotSupport))
            d['intron_annotation_support'] = ','.join(map(str, s.IntronAnnotSupport))
            metrics['Splice Support']['unknown_likely_coding'].append(s.IntronRnaSupportPercent)
            metrics['Exon Support']['unknown_likely_coding'].append(s.ExonRnaSupportPercent)


def validate_pacbio_splices(deduplicated_strand_resolved_consensus, db_path, tx_dict, metrics, require_pacbio_support):
    """
    Tag transcripts as having PacBio support.
    If users passed the --require-pacbio-support, remove any transcript which does not have support.

    TODO: consider doing fuzzy matching due to noisy nature of PacBio reads.
    """
    pb_intervals = tools.sqlInterface.load_pb_intron_intervals(db_path)
    pb_resolved_consensus = []
    for tx_id, d in deduplicated_strand_resolved_consensus:
        tx = tx_dict[tx_id]
        # remove strand information from the existing intervals
        intervals = frozenset([tools.intervals.ChromosomeInterval(i.chromosome, i.start, i.stop, '.')
                               for i in tx.intron_intervals])
        if intervals in pb_intervals:
            d['pacbio_isoform_supported'] = True
            metrics['IsoSeq Transcript Valdiation'][True] += 1
            pb_resolved_consensus.append([tx_id, d])
        elif require_pacbio_support is False:
            d['pacbio_isoform_supported'] = False
            metrics['IsoSeq Transcript Valdiation'][False] += 1
            pb_resolved_consensus.append([tx_id, d])
        # if require_pacbio_support is True, then we don't save this transcript
    return pb_resolved_consensus


def is_failed_df(df):
    """Failed genes have no passing transcripts. Handles series"""
    try:
        return not df.TranscriptClass.str.contains('passing').any()
    except AttributeError:
        return df.TranscriptClass == 'failing'


def rescue_failed_gene(gene_df, gene_id, metrics, hints_db_has_rnaseq):
    """Rescues a failed gene by picking the one transcript with the highest coverage"""
    gene_df = gene_df.sort_values('AlnCoverage', ascending=False)
    best_rows = find_best_score(gene_df, column='AlnCoverage')
    return incorporate_tx(best_rows, gene_id, metrics, hints_db_has_rnaseq, failed_gene=True)


def incorporate_tx(best_rows, gene_id, metrics, hints_db_has_rnaseq, failed_gene):
    """incorporate a transcript into the consensus set, storing metrics."""
    best_series = best_rows.iloc[0]
    transcript_modes = evaluate_ties(best_rows)
    # construct the tags for this transcript
    d = {'source_transcript': best_series.name,
         'source_gene': gene_id,
         'score': round(best_series.AlnGoodness, 3),
         'failed_gene': failed_gene,
         'transcript_modes': transcript_modes,
         'gene_biotype': best_series.GeneBiotype,
         'transcript_class': best_series.TranscriptClass,
         'transcript_biotype': best_series.TranscriptBiotype,
         'exon_annotation_support': ','.join(map(str, best_series.ExonAnnotSupport)),
         'intron_annotation_support': ','.join(map(str, best_series.IntronAnnotSupport))}
    if best_series.TranscriptBiotype == 'protein_coding':
        d['cds_annotation_support'] = ','.join(map(str, best_series.CdsAnnotSupport))
    if hints_db_has_rnaseq is True:
        d['exon_rna_support'] = ','.join(map(str, best_series.ExonRnaSupport))
        d['intron_rna_support'] = ','.join(map(str, best_series.IntronRnaSupport))
    if best_series.Paralogy > 1:
        assert best_series.ParalogStatus is not None
        d['paralogy'] = best_series.Paralogy
        d['paralog_status'] = best_series.ParalogStatus
    if 'GeneAlternateContigs' in best_series and best_series.GeneAlternateContigs is not None:
        d['gene_alterate_contigs'] = best_series.GeneAlternateContigs
    if best_series.GeneName is not None:
        d['source_gene_common_name'] = best_series.GeneName

    # add information to the overall metrics
    if best_series.TranscriptBiotype == 'protein_coding':
        metrics['Transcript Modes'][transcript_modes] += 1
    metrics['Coverage'][best_series.TranscriptBiotype].append(best_series.AlnCoverage)
    metrics['Identity'][best_series.TranscriptBiotype].append(best_series.AlnIdentity)
    metrics['Splice Support'][best_series.TranscriptBiotype].append(best_series.IntronRnaSupportPercent)
    metrics['Exon Support'][best_series.TranscriptBiotype].append(best_series.ExonRnaSupportPercent)
    metrics['Splice Annotation Support'][best_series.TranscriptBiotype].append(best_series.IntronAnnotSupportPercent)
    metrics['Exon Annotation Support'][best_series.TranscriptBiotype].append(best_series.ExonAnnotSupportPercent)
    metrics['Original Introns'][best_series.TranscriptBiotype].append(best_series.OriginalIntronsPercent)
    return best_series.AlignmentId, d


def find_best_score(tx_df, column='TranscriptScore'):
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


def evaluate_ties(best_rows):
    """Find out how many transcript modes agreed on this"""
    return ','.join(sorted(set([tools.nameConversions.alignment_type(x) for x in best_rows.AlignmentId])))


def slice_df(df, ix):
    """
    Slices a DataFrame by an index, handling the case where the index is missing. CHandles the case where a single row
    is returned, thus making it a series.
    """
    try:
        r = df.xs(ix)
        if isinstance(r, pd.core.series.Series):
            return pd.DataFrame([r])
        else:
            return r
    except KeyError:
        return pd.DataFrame()


def find_novel_splices(gene_consensus_dict, denovo_gene_df, tx_dict, gene_id, common_name_map, metrics, failed_gene,
                       gene_biotype, denovo_num_introns):
    """
    Finds novel splice junctions in CGP/PB transcripts. A novel splice junction is defined as a splice which
    homGeneMapping did not map over and which is supported by RNA-seq.
    """
    # extract all splices we have already seen for this gene
    existing_splices = set()
    for consensus_tx in gene_consensus_dict:
        existing_splices.update(tx_dict[consensus_tx].intron_intervals)

    denovo_tx_dict = {}
    for aln_id, s in denovo_gene_df.iterrows():
        tx_mode = tools.nameConversions.alignment_type(aln_id)
        denovo_tx_obj = tx_dict[aln_id]
        if len(denovo_tx_obj.intron_intervals) < denovo_num_introns:
            continue
        new_supported_splices = set()
        for intron, rna in zip(*[denovo_tx_obj.intron_intervals, s.IntronRnaSupport]):
            if rna > 0 and intron not in existing_splices:
                new_supported_splices.add(intron)
        if len(new_supported_splices) == 0:
            continue  # nothing novel here
        # if any splices are both not supported by annotation and supported by RNA, call this as novel
        if any(annot == 0 and i in new_supported_splices for i, annot in zip(*[denovo_tx_obj.intron_intervals,
                                                                               s.IntronAnnotSupport])):
            metrics['denovo'][tx_mode]['Novel isoforms'] += 1
            metrics['Transcript Modes'][tx_mode] += 1
            tx_class = 'putative_novel_isoform'
        # if any splices are new, and supported by RNA-seq call this poor alignment
        else:
            metrics['Transcript Modes'][tx_mode] += 1
            tx_class = 'poor_alignment'
        denovo_tx_dict[aln_id] = {'transcript_class': tx_class, 'source_gene': gene_id, 'failed_gene': failed_gene,
                                  'transcript_biotype': 'unknown_likely_coding', 'gene_biotype': gene_biotype,
                                  'intron_rna_support': ','.join(map(str, s.IntronRnaSupport)),
                                  'exon_rna_support': ','.join(map(str, s.ExonRnaSupport)),
                                  'transcript_modes': tx_mode,
                                  'exon_annotation_support': ','.join(map(str, s.ExonAnnotSupport)),
                                  'intron_annotation_support': ','.join(map(str, s.IntronAnnotSupport)),
                                  'cds_annotation_support': ','.join(map(str, s.CdsAnnotSupport))}
        common_name = common_name_map[gene_id]
        if common_name != gene_id:
            denovo_tx_dict[aln_id]['source_gene_common_name'] = common_name
        metrics['Splice Support']['unknown_likely_coding'].append(s.IntronRnaSupportPercent)
        metrics['Exon Support']['unknown_likely_coding'].append(s.ExonRnaSupportPercent)
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
                               key=lambda (tx, s): -s)
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
        # don't count novel transcripts towards completeness
        if tools.nameConversions.aln_id_is_cgp(aln_id) or tools.nameConversions.aln_id_is_pb(aln_id):
            continue
        genes[c['gene_biotype']].add(c['source_gene'])
        txs[c['transcript_biotype']] += 1
    genes = {biotype: len(gene_list) for biotype, gene_list in genes.iteritems()}
    metrics['Completeness'] = {'Gene': genes, 'Transcript': txs}


def calculate_improvement_metrics(final_consensus, scored_df, tm_eval_df, metrics):
    """For coding transcripts, how much did we improve the metrics?"""
    df = scored_df.reset_index()
    df = df.set_index('AlignmentId')
    tm_df = tm_eval_df.reset_index()
    tm_df = tm_df.set_index('AlignmentId')
    metrics['Evaluation Improvement'] = []
    for aln_id, c in final_consensus:
        if c['TranscriptBiotype'] != 'protein_coding':
            continue
        s = df.ix[aln_id]
        tm_s = tm_df.ix[tools.nameConversions.remove_augustus_alignment_number(aln_id)]
        metrics['Evaluation Improvement'].append([tm_s.TransMapOriginalIntronsPercent,
                                                  tm_s.IntronAnnotSupportPercent,
                                                  tm_s.IntronRnaSupportPercent,
                                                  s.OriginalIntronsPercent,
                                                  s.IntronAnnotSupportPercent,
                                                  s.IntronRnaSupportPercent])


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
            tx_obj.score = attrs.get('score', 0)
            tx_count += 1
            source_gene = attrs.get('source_gene', tx_obj.name2)
            if source_gene not in genes_seen:
                genes_seen.add(source_gene)
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
        if 'score' in attrs:
            score = attrs['score']
            del attrs['score']
        else:
            score = '.'
        if 'source_gene_common_name' in attrs:
            attrs['Name'] = attrs['source_gene_common_name']
        # don't include the support vectors in the string, they will be placed in their respective places
        attrs_str = ['='.join([key, str(val)]) for key, val in sorted(attrs.iteritems()) if 'support' not in key]
        return score, ';'.join(attrs_str)

    def find_feature_support(attrs, feature, i):
        """Extracts the boolean value from the comma delimited string"""
        vals = map(bool, attrs[feature].split(','))
        return vals[i]

    def generate_gene_record(chrom, tx_objs, gene_id, attrs):
        """calculates the gene interval for this list of tx"""
        intervals = set()
        for tx in tx_objs:
            intervals.update(tx.exon_intervals)
        intervals = sorted(intervals)
        strand = tx_objs[0].strand
        # subset the attrs to gene fields
        useful_keys = ['source_gene_common_name', 'source_gene', 'gene_biotype', 'failed_gene', 'transcript_modes',
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
        for line in generate_intron_exon_records(chrom, tx_obj, tx_id, attrs):
            yield line
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
            yield [chrom, 'CAT', 'exon', exon.start + 1, exon.stop + 1, score, exon.strand, '.', attrs_field]
            cds_interval = exon.intersection(tx_obj.coding_interval)
            if cds_interval is not None:
                attrs['reference_support'] = find_feature_support(attrs, 'cds_annotation_support', cds_i)
                score, attrs_field = convert_attrs(attrs, 'CDS:{}:{}'.format(tx_id, cds_i))
                cds_i += 1
                yield [chrom, 'CAT', 'CDS', cds_interval.start + 1, cds_interval.stop + 1, score, exon.strand,
                       convert_frame(exon_frame), attrs_field]

        # intron records
        for i, intron in enumerate(tx_obj.intron_intervals):
            attrs['rna_support'] = find_feature_support(attrs, 'intron_rna_support', i)
            attrs['reference_support'] = find_feature_support(attrs, 'intron_annotation_support', i)
            score, attrs_field = convert_attrs(attrs, 'exon:{}:{}'.format(tx_id, i))
            yield [chrom, 'CAT', 'intron', intron.start + 1, intron.stop + 1, score, intron.strand, '.', attrs_field]

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
