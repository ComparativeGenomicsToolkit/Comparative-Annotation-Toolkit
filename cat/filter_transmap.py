"""
Resolves paralogs in transMap output based on MLE estimate of the distribution of alignment identities in the transMap
process.

Resolves genes that have transcripts split across chromosomes based on a consensus finding process. This process
combines information from both synteny and phylogenetic distance to determine which contig is likely the parental
contig. This should be used carefully on genomes that are highly fragmented.

"""
import logging
import collections
import numpy as np
import pandas as pd
from scipy.stats import norm
import tools.nameConversions
import tools.sqlInterface
import tools.transcripts
import tools.mathOps
import tools.fileOps
import tools.intervals

pd.options.mode.chained_assignment = None
logger = logging.getLogger(__name__)


def filter_transmap(filter_tm_args, out_target):
    """
    Entry point for transMap filtering.
    :param filter_tm_args: argparse Namespace produced by FilterTransMap.get_args()
    :param out_target: luigi.LocalTarget where the results will be written
    """
    # load database tables
    ref_df = tools.sqlInterface.load_annotation(filter_tm_args.ref_db_path)
    aln_eval_df = tools.sqlInterface.load_alignment_evaluation(filter_tm_args.db_path)
    tx_dict = tools.transcripts.get_gene_pred_dict(filter_tm_args.tm_gp)

    # fit lognormal distributions to assign categories
    updated_aln_eval_df, fit_df = fit_distributions(aln_eval_df, ref_df, filter_tm_args.genome)

    # store metrics on filtering for plots. This will be written to disk for final plots
    metrics = {}

    # resolve paralogs
    paralog_metrics, paralog_resolved_df = resolve_paralogs(updated_aln_eval_df, filter_tm_args.genome)
    metrics['Paralogy'] = paralog_metrics

    # resolve split genes. If user requests the bad ones to be removed, do so. Make an attempt to rescue.
    split_gene_metrics, resolved_df = resolve_split_genes(paralog_resolved_df, tx_dict,
                                                          filter_tm_args.resolve_split_genes,
                                                          filter_tm_args.genome)
    metrics['Split Genes'] = split_gene_metrics

    # filter the paralogous mappings left
    paralog_filtered_df = resolved_df[resolved_df.UpdatedParalogStatus != 'ToRemove']
    paralog_filtered_df = paralog_filtered_df.drop('ParalogStatus', axis=1)
    paralog_filtered_df = paralog_filtered_df.rename(columns={'UpdatedParalogStatus': 'ParalogStatus'})
    assert len(paralog_filtered_df) == len(set(paralog_filtered_df.TranscriptId))

    # keep some metrics
    for biotype, biotype_df in resolved_df.groupby('TranscriptBiotype'):
        metrics['Paralogy'][biotype]['Alignments discarded'] = len(biotype_df[biotype_df.UpdatedParalogStatus == 'ToRemove'])
        tot = paralog_metrics[biotype]['Synteny heuristic'] + paralog_metrics[biotype]['Model prediction'] + \
              paralog_metrics[biotype]['Arbitrarily resolved']
        logger.info('Discarded {:,} alignments for {} on {} after paralog resolution. '
                    '{:,} transcripts remain.'.format(paralog_metrics[biotype]['Alignments discarded'],
                                                      biotype, filter_tm_args.genome, tot))

    # write out the filtered transMap results
    with out_target.open('w') as outf:
        for aln_id in paralog_filtered_df.AlignmentId:
            tx = tx_dict[aln_id]
            outf.write('\t'.join(tx.get_gene_pred()) + '\n')

    # produced update df
    updated_df = create_new_table(paralog_filtered_df)

    return metrics, updated_df, fit_df


def fit_distributions(aln_eval_df, ref_df, genome):
    """
    Fits a normal distribution to the -log(100 - identity) where identity != 100. Uses the MLE estimate to determine
    a cutoff of identity specific to this genetic distance.
    """
    def transform_data(idents):
        """transforms identity data to -log(100 - ident) where ident != 100"""
        return -np.log(100 - idents[idents != 100])

    def find_cutoff(biotype_unique, num_sigma=1):
        """Locates the MLE identity cutoff"""
        unique_mu, unique_sigma = norm.fit(transform_data(biotype_unique.TransMapIdentity))
        return 100 - np.exp(-(unique_mu - (num_sigma * unique_sigma)))

    biotype_df = pd.merge(aln_eval_df, ref_df, on='TranscriptId')
    r = []  # will hold the labels
    identity_cutoffs = []

    for biotype, df in biotype_df.groupby('TranscriptBiotype'):
        biotype_unique = df[df.Paralogy == 1]
        biotype_not_unique = df[df.Paralogy != 1]
        if len(biotype_not_unique) == 0:
            # No paralogous mappings implies all passing
            logger.info('No paralogous mappings for {} on {}.'.format(biotype, genome))
            r.extend([[aln_id, 'passing'] for aln_id in df.AlignmentId])
            identity_cutoffs.append([biotype, None])
        elif len(biotype_unique) == 0:
            # Only paralogous mappings implies all failing
            logger.info('Only paralogous mappings for {} on {}.'.format(biotype, genome))
            r.extend([[aln_id, 'failing'] for aln_id in df.AlignmentId])
            identity_cutoffs.append([biotype, None])
        else:
            logger.info('{:,} paralogous mappings and {:,} 1-1 orthologous mappings'
                        ' for {} on {}.'.format(len(biotype_not_unique), len(biotype_unique), biotype, genome))
            cutoff = find_cutoff(biotype_unique)
            identity_cutoffs.append([biotype, cutoff])
            if np.isnan(cutoff) or np.isinf(cutoff):
                logger.warning('Unable to establish a identity boundary for {} on {}. '
                               'All transcripts marked passing.'.format(biotype, genome))
                labels = ['passing'] * len(df)
            else:
                labels = ['passing' if ident >= cutoff else 'failing' for ident in df.TransMapIdentity]
                num_pass = labels.count('passing')
                num_fail = labels.count('failing')
                logger.info('Established a {:.2%} identity boundary for {} on {} resulting in '
                            '{:,} passing and {:,} failing alignments.'.format(cutoff / 100, biotype, genome,
                                                                               num_pass, num_fail))
            r.extend([[aln_id, label] for aln_id, label in zip(*[df.AlignmentId, labels])])

    # turn r into a dataframe
    r_df = pd.DataFrame(r)
    r_df.columns = ['AlignmentId', 'TranscriptClass']
    ident_df = pd.DataFrame(identity_cutoffs)
    ident_df.columns = ['TranscriptBiotype', 'IdentityCutoff']
    return pd.merge(biotype_df, r_df, on='AlignmentId'), ident_df


def resolve_paralogs(updated_aln_eval_df, genome):
    """
    Resolves paralogs based on likelihood to come from the two lognormal distributions.

    1. If only one paralog is more likely under the ortholog model, discard the others.
    2. If more than one paralog are more likely under the ortholog model, or we were unable to fit a model, resolve
       based on the synteny score. See calculate_synteny_score

    :param updated_aln_eval_df: DataFrame produced by fit_distributions()
    :return: tuple of (metrics_dict, filtered DataFrame)
    """
    def apply_label(s):
        return s.TranscriptClass if s.ParalogStatus != 'NotConfident' else 'failing'

    updated_aln_eval_df['Score'] = updated_aln_eval_df.apply(calculate_synteny_score, axis=1)
    updated_aln_eval_df = updated_aln_eval_df.sort_values(by='Score', ascending=False)

    paralog_status = []  # stores the results for a new column
    paralog_metrics = {biotype: {'Alignments discarded': 0, 'Model prediction': 0,
                                 'Synteny heuristic': 0, 'Arbitrarily resolved': 0}
                       for biotype in set(updated_aln_eval_df.TranscriptBiotype)}

    for tx, df in updated_aln_eval_df.groupby('TranscriptId'):
        if len(df) == 1:  # no paralogs
            paralog_status.append([df.AlignmentId.iloc[0], 'Unique'])
            continue
        biotype = df.TranscriptBiotype.iloc[0]
        passing = df[df.TranscriptClass == 'passing']
        if len(passing) == 1:  # we can pick one passing member
            paralog_metrics[biotype]['Model prediction'] += 1
            paralog_status.append([df.AlignmentId.iloc[0], 'ModelPrediction'])
        else:
            highest_score_df = df[df.Score == df.iloc[0].Score]
            if len(highest_score_df) == 1:
                paralog_metrics[biotype]['Synteny heuristic'] += 1
                paralog_status.append([highest_score_df.AlignmentId.iloc[0], 'SyntenyHeuristic'])
            else:
                paralog_metrics[biotype]['Arbitrarily resolved'] += 1
                paralog_status.append([highest_score_df.AlignmentId.iloc[0], 'NotConfident'])

    # keep track of the ones we flagged for split genes
    kept_ids = set(zip(*paralog_status)[0])
    all_ids = set(updated_aln_eval_df.AlignmentId)
    paralog_status.extend([[aln_id, 'ToRemove'] for aln_id in all_ids - kept_ids])

    status_df = pd.DataFrame(paralog_status)
    status_df.columns = ['AlignmentId', 'ParalogStatus']
    merged = pd.merge(status_df, updated_aln_eval_df, on='AlignmentId')
    merged['TranscriptClass'] = merged.apply(apply_label, axis=1)
    return paralog_metrics, merged


def resolve_split_genes(paralog_resolved_df, tx_dict, remove_split_genes, genome):
    """
    Resolves cases where transMap mapped a gene to different chromosomes. This is a useful feature to turn on
    if you have a high quality assembly, but may be problematic for highly fragmented assemblies.

    For each gene, find all transcript clusters. Pick the cluster with the highest average synteny score. Keep track
    of whether we split over contigs and/or over different regions of the same contig.

    :param paralog_resolved_df: DataFrame produced by resolve_paralogs()
    :param tx_dict: Dictionary mapping alignment IDs to GenePredTranscript objects.
    :param remove_split_genes: Boolean. Do we remove the transcripts on the less likely contig after resolution?
    :param genome: Genome in question. For logging.
    :return: tuple of (metrics_dict, updated dataframe)
    """
    def cluster_gene(tx_objs, paralog_statuses):
        """
        Cluster together all transcripts of a given gene
        """
        # split up the transcripts by chromosome to make merging easier
        intervals_by_chrom = collections.defaultdict(list)
        for tx_obj in tx_objs:
            intervals_by_chrom[tx_obj.chromosome].append(tx_obj.interval)

        # cluster all intervals for a chromosome
        merged_intervals = []
        for chrom, intervals in intervals_by_chrom.iteritems():
            merged_intervals.extend(tools.intervals.gap_merge_intervals(intervals, 0))

        # assign each original tx_obj to a cluster
        clusters = collections.defaultdict(list)
        for tx_obj, paralog_status in zip(*[tx_objs, paralog_statuses]):
            for i, interval in enumerate(merged_intervals):
                if interval.overlap(tx_obj.interval):
                    clusters[i].append([tx_obj.name, paralog_status])
                    break

        # remove any clusters that contain only paralogous transcripts
        filtered_clusters = collections.defaultdict(list)
        for i, cluster in clusters.iteritems():
            for aln_id, paralog_status in cluster:
                if paralog_status != 'ToRemove':
                    filtered_clusters[i].append(aln_id)
        return clusters, filtered_clusters

    def find_best_cluster(clusters, rec):
        """
        If we have more than one cluster for a gene, try to resolve it. Report on the results.
        """
        scores = {}
        for cluster_id, cluster_tx_ids in clusters.iteritems():
            scores[cluster_id] = rec[rec.AlignmentId.isin(cluster_tx_ids)].Score.mean()
        sorted_scores = sorted(scores.iteritems(), key=lambda (cid, s): -s)
        return sorted_scores[0][0]

    def is_split_chrom_gene(tx_objs):
        """
        If we are trying to resolve a split gene, is it split over contigs or just location?
        """
        return len({tx_obj.chromosome for tx_obj in tx_objs}) != 1

    def is_split_same_chrom_gene(gene_tx_obj_dict, clusters):
        """
        If we are trying to resolve a split gene, is it split over the same contig?
        """
        cluster_chroms = set()
        for cluster_id, cluster_tx_ids in clusters.iteritems():
            chrom = gene_tx_obj_dict[cluster_tx_ids[0]].chromosome
            if chrom in cluster_chroms:
                return True
            cluster_chroms.add(chrom)
        return False

    split_gene_metrics = {'Number of contig split genes': 0,
                          'Number of intra-contig split genes': 0}
    alignment_ids_to_remove = set()
    split_status = []
    for gene, rec in paralog_resolved_df.groupby('GeneId'):
        gene_tx_obj_dict = collections.OrderedDict([[tx_id, tx_dict[tx_id]] for tx_id in rec.AlignmentId])
        paralog_statuses = collections.OrderedDict([[tx_id, s] for tx_id, s in zip(*[rec.AlignmentId, rec.ParalogStatus])])
        clusters, filtered_clusters = cluster_gene(gene_tx_obj_dict.values(), paralog_statuses.values())

        if len(filtered_clusters) == 1:  # no split contigs here
            split_status.extend([[aln_id, None, False, s] for aln_id, s in zip(*[rec.AlignmentId, rec.ParalogStatus])])
            continue

        # find the best cluster based on score
        best_cluster = find_best_cluster(filtered_clusters, rec)

        # extract IDs that do not belong to the best cluster, only if they are not already filtered
        for cluster, cluster_ids in clusters.iteritems():
            if cluster == best_cluster:
                continue
            for aln_id, paralog_status in cluster_ids:
                if paralog_status != 'ToRemove':
                    alignment_ids_to_remove.add(aln_id)

        # record some metrics
        if is_split_chrom_gene(gene_tx_obj_dict.values()):
            split_gene_metrics['Number of contig split genes'] += 1
            if is_split_same_chrom_gene(gene_tx_obj_dict, filtered_clusters) is True:
                split_gene_metrics['Number of intra-contig split genes'] += 1
        else:
            split_gene_metrics['Number of intra-contig split genes'] += 1

        # begin iterating over clusters, assigning alternative chromosomes
        for cluster, cluster_ids in clusters.iteritems():
            chrom = gene_tx_obj_dict[cluster_ids[0][0]].chromosome
            other_chroms = {tx_obj.chromosome for tx_obj in gene_tx_obj_dict.itervalues() if tx_obj.chromosome != chrom}
            other_chroms = ','.join(other_chroms)
            # begin paralog rescue, dealing with the rare case where we have chaining problems
            counts = collections.Counter([tools.nameConversions.strip_alignment_numbers(aln_id)
                                          for aln_id, s in cluster_ids])
            for aln_id, paralog_status in cluster_ids:
                tx_id = tools.nameConversions.strip_alignment_numbers(aln_id)
                if paralog_status == 'ToRemove' and cluster == best_cluster and counts[tx_id] == 1:
                    paralog_status = 'Rescued'
                split_status.append([aln_id, other_chroms, True, paralog_status])

    split_df = pd.DataFrame(split_status)
    split_df.columns = ['AlignmentId', 'GeneAlternateContigs', 'SplitGene', 'UpdatedParalogStatus']
    split_gene_metrics['Number of transcripts rescued'] = len(split_df[split_df.UpdatedParalogStatus == 'Rescued'])
    merged_df = pd.merge(paralog_resolved_df, split_df, on='AlignmentId')

    if remove_split_genes:
        split_gene_metrics['Number of transcripts removed'] = len(alignment_ids_to_remove)
        logger.info('{:,} genes for {} have transcripts split across contigs. '
                    '{:,} transcripts removed.'.format(split_gene_metrics['Number of contig split genes'],
                                                       genome,
                                                       split_gene_metrics['Number of transcripts removed']))
        filtered_df = merged_df[~merged_df['AlignmentId'].isin(alignment_ids_to_remove)]
        return split_gene_metrics, filtered_df
    else:
        logger.info('{:,} genes for {} have transcripts '
                    'split across contigs.'.format(split_gene_metrics['Number of contig split genes'], genome))
        return split_gene_metrics, merged_df


def create_new_table(paralog_filtered_df):
    """
    Clean up this dataframe to write to SQL. I am not using the long form here because it is too hard to reload when
    the columns are not numeric.
    :param paralog_filtered_df: output from resolve_split_genes()
    :return: dataframe to be written to sql
    """
    cols = ['GeneId', 'TranscriptId', 'AlignmentId', 'TranscriptClass', 'ParalogStatus', 'GeneAlternateContigs',
            'SplitGene']
    df = paralog_filtered_df[cols]
    return df.set_index(['TranscriptId', 'AlignmentId'])


def calculate_synteny_score(s):
    """
    Function to score an alignment. Scoring method:
    0.2 * coverage + 0.3 * identity + 0.5 * synteny
    :param s: pandas Series
    :return: float between 0 and 100
    """
    r = 0.2 * s.TransMapCoverage + \
        0.3 * s.TransMapIdentity + \
        0.5 * (1.0 * s.Synteny / 10)
    assert 0 <= r <= 100
    return r
