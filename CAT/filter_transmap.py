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
    paralog_metrics, paralog_filtered_df = resolve_paralogs(updated_aln_eval_df, filter_tm_args.genome)
    metrics['Paralogy'] = paralog_metrics

    # resolve split genes, if user requested
    if filter_tm_args.resolve_split_genes is True:
        split_gene_metrics, paralog_filtered_df = resolve_split_genes(paralog_filtered_df, tx_dict)
        metrics['Split Genes'] = split_gene_metrics

    # write out the filtered transMap results
    with out_target.open('w') as outf:
        for aln_id in paralog_filtered_df.AlignmentId:
            tx = tx_dict[aln_id]
            outf.write('\t'.join(tx.get_gene_pred()) + '\n')

    # produced update df
    updated_df = create_new_table(paralog_filtered_df, filter_tm_args.resolve_split_genes)

    return metrics, updated_df, fit_df


def fit_distributions(aln_eval_df, ref_df, genome):
    """
    Fits a normal distribution to the -log(1 - identity) where identity != 1. Uses the MLE estimate to determine
    a cutoff of identity specific to this genetic distance.
    """
    def transform_data(idents):
        """transforms identity data to -log(1 - ident) where ident != 1"""
        return -np.log(1 - idents[idents != 1])

    def find_cutoff(biotype_unique, num_sigma=1):
        """Locates the MLE identity cutoff"""
        unique_mu, unique_sigma = norm.fit(transform_data(biotype_unique.TransMapIdentity))
        return 1 - np.exp(-(unique_mu - (num_sigma * unique_sigma)))

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
            labels = ['passing' if ident >= cutoff else 'failing' for ident in df.TransMapIdentity]
            num_pass = labels.count('passing')
            num_fail = labels.count('failing')
            logger.info('Established a {:.2%} identity boundary for {} on {} resulting in '
                        '{:,} passing and {:,} failing alignments.'.format(cutoff, biotype, genome, num_pass, num_fail))
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
            paralog_status.append([df.AlignmentId.iloc[0], None])
            continue
        biotype = df.TranscriptBiotype.iloc[0]
        passing = df[df.TranscriptClass == 'passing']
        if len(passing) == 1:  # we can pick one passing member
            paralog_metrics[biotype]['Model prediction'] += 1
            paralog_metrics[biotype]['Alignments discarded'] += len(df) - 1
            paralog_status.append([df.AlignmentId.iloc[0], 'Confident'])
        else:
            highest_score_df = df[df.Score == df.iloc[0].Score]
            if len(highest_score_df) == 1:
                paralog_metrics[biotype]['Synteny heuristic'] += 1
                paralog_metrics[biotype]['Alignments discarded'] += len(df) - 1
                paralog_status.append([highest_score_df.AlignmentId.iloc[0], 'Confident'])
            else:
                paralog_metrics[biotype]['Arbitrarily resolved'] += 1
                paralog_metrics[biotype]['Alignments discarded'] += len(df) - 1
                paralog_status.append([highest_score_df.AlignmentId.iloc[0], 'NotConfident'])

    for biotype in set(updated_aln_eval_df.TranscriptBiotype):
        tot = paralog_metrics[biotype]['Synteny heuristic'] + paralog_metrics[biotype]['Model prediction'] + \
              paralog_metrics[biotype]['Arbitrarily resolved']
        logger.info('Discarded {:,} alignments for {} on {} after paralog resolution. '
                    '{:,} transcripts remain.'.format(paralog_metrics[biotype]['Alignments discarded'],
                                                      biotype, genome, tot))
    status_df = pd.DataFrame(paralog_status)
    status_df.columns = ['AlignmentId', 'ParalogStatus']
    merged = pd.merge(status_df, updated_aln_eval_df, on='AlignmentId')  # this filters out paralogous alignments
    merged['TranscriptClass'] = merged.apply(apply_label, axis=1)
    return paralog_metrics, merged


def resolve_split_genes(paralog_filtered_df, tx_dict):
    """
    Resolves cases where transMap mapped a gene to different chromosomes. This is a useful feature to turn on
    if you have a high quality assembly, but may be problematic for highly fragmented assemblies.

    For each gene, if transcripts on that gene are on multiple sequences, a consensus finding process is performed.
    For each chromosome a gene maps to, calculate the same synteny metric in paralog resolution, and then find
    the chromosome with the highest average metric.

    First look only at transcripts which match the gene biotype.

    :param paralog_filtered_df: DataFrame produced by resolve_paralogs()
    :param tx_dict: Dictionary mapping alignment IDs to GenePredTranscript objects.
    :return: tuple of (metrics_dict, updated dataframe)
    """
    def extract_gene_biotype(rec):
        """there should only be one gene biotype"""
        gene_biotypes = list(set(rec.GeneBiotype))
        assert len(gene_biotypes) == 1
        return gene_biotypes[0]

    def chrom_metric(chroms, rec):
        """Finds the average score for each chromosome scored"""
        r = {}
        tot = len(rec)
        for chrom, vals in chroms.iteritems():
            r[chrom] = (1.0 * len(vals) / tot) * np.mean(zip(*vals)[1])
        return r

    def find_chromosomes(rec, gene_biotype=None):
        """create a mapping of genome sequence names to associated alignment ids and synteny scores"""
        chroms = collections.defaultdict(list)
        for _, s in rec.iterrows():
            if gene_biotype is not None and s.TranscriptBiotype != gene_biotype:
                continue
            chroms[tx_dict[s.AlignmentId].chromosome].append([s.AlignmentId, calculate_synteny_score(s)])
        return chroms

    def find_best_chroms(chroms):
        """finds the best chrom in the chroms dict based on the highest score score"""
        chroms = chrom_metric(chroms, rec)
        s = sorted(chroms.iteritems(), key=lambda (chrom, val): val)
        best_val = s[0][1]
        return [chrom for chrom, val in s if val == best_val]

    def find_names_to_remove(chrom_or_chroms, chroms):
        """return the list of names to remove as well as the count based on chrom_or_chroms"""
        if isinstance(chrom_or_chroms, str):
            chrom_list = [chrom_or_chroms]
        else:
            chrom_list = chrom_or_chroms
        aln_ids = set()
        for chrom in chrom_list:
            aln_ids.update([aln_id for aln_id, score in chroms[chrom]])
        return len(aln_ids), aln_ids

    split_gene_metrics = {'Number of split genes': 0, 'Number of transcripts removed': 0}
    alignment_ids_to_remove = set()
    split_status = []
    for gene, rec in paralog_filtered_df.groupby('GeneId'):
        gene_biotype = extract_gene_biotype(rec)
        if gene_biotype == 'protein_coding':
            chroms = find_chromosomes(rec, gene_biotype)
            if len(chroms) == 0:  # there may be an edge case where a protein coding gene has no coding transcripts
                chroms = find_chromosomes(rec)
        else:
            chroms = find_chromosomes(rec)
        if len(chroms) == 1:  # no ambiguous chromosomes here
            split_status.extend([[tx_id, None] for tx_id in rec.TranscriptId])
            continue
        best_chroms = find_best_chroms(chroms)
        n, i = find_names_to_remove(best_chroms[0], chroms)
        alignment_ids_to_remove.update(i)
        split_gene_metrics['Number of split genes'] += 1
        split_gene_metrics['Number of transcripts removed'] += len(i)
        other_contigs = ','.join(set(chroms.keys()) - {best_chroms[0]})
        split_status.extend([[tx_id, other_contigs] for tx_id in rec.TranscriptId])

    split_df = pd.DataFrame(split_status)
    split_df.columns = ['TranscriptId', 'GeneAlternateContigs']
    merged_df = pd.merge(paralog_filtered_df, split_df, on='TranscriptId')
    final_df = merged_df[~merged_df['AlignmentId'].isin(alignment_ids_to_remove)]
    return split_gene_metrics, final_df


def create_new_table(paralog_filtered_df, resolve_split_genes_flag):
    """
    Clean up this dataframe to write to SQL. I am not using the long form here because it is too hard to reload when
    the columns are not numeric.
    :param paralog_filtered_df: output from either resolve_split_genes() if flag set else resolve_paralogs()
    :return: dataframe to be written to sql
    """
    cols = ['GeneId', 'TranscriptId', 'AlignmentId', 'TranscriptClass', 'ParalogStatus']
    if resolve_split_genes_flag is True:
        cols.append('GeneAlternateContigs')
    df = paralog_filtered_df[cols]
    return df.set_index(['TranscriptId', 'AlignmentId'])


def calculate_synteny_score(s):
    """
    Function to score an alignment. Scoring method:
    0.2 * coverage + 0.3 * identity + 0.5 * synteny
    :param s: pandas Series
    :return: float between 0 and 1
    """
    r = 0.002 * s.TransMapCoverage + \
        0.003 * s.TransMapIdentity + \
        0.5 * (1.0 * s.Synteny / 6)
    assert 0 <= r <= 1
    return r
