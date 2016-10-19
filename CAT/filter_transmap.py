"""
Resolves paralogs in transMap output based on phylogenetic distance inferred by FastTree.
Resolves genes that have transcripts split across chromosomes based on a consensus finding process. This process
combines information from both synteny and phylogenetic distance to determine which contig is likely the parental
contig. This should be used carefully on genomes that are highly fragmented.
"""
import collections
import numpy as np
import pandas as pd
from pomegranate import NormalDistribution, GeneralMixtureModel
from scipy.stats import norm
import tools.nameConversions
import tools.sqlInterface
import tools.transcripts
import tools.mathOps

pd.options.mode.chained_assignment = None


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
    updated_aln_eval_df, fit_df = fit_distributions(aln_eval_df, ref_df)

    # store metrics on filtering for plots. This will be written to disk for final plots
    metrics = {}

    # resolve paralogs
    paralog_metrics, paralog_filtered_df = resolve_paralogs(updated_aln_eval_df)
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
    updated_df = create_new_table(paralog_filtered_df)

    return metrics, updated_df, fit_df


def fit_distributions(aln_eval_df, ref_df):
    """
    Fits a mixture model of 2 lognormals to the aln_eval_df identity based on whether the alignments are 1-1 or 1-many
    This is done for every biotype. In the end, a new column is added to the aln_eval_df that says whether a alignment
    is Failing (more likely under paralog model) or Passing (more likely under paralog model)
    """
    def transform_data(idents):
        """transforms identity data to -log(1 - ident) where ident != 1"""
        return -np.log(1 - idents[idents != 1])

    def reverse_transform(x):
        """transforms an identity in reverse back to % id"""
        return 1 - np.exp(-x)

    def calculate_starting_values(idents):
        """calculates normal distribution parameters on the quantity log(1 - identity) where identity != 1"""
        return norm.fit(idents)

    def find_xvals(mu, sigma):
        """infers the likely useful xvalues by using the interval function to find where 99% of the data are"""
        start, stop = norm.interval(0.99, mu, sigma)
        return np.floor(start), np.ceil(stop)

    def find_cutoff(d):
        """Locates the identity cutoff point"""
        unique_start, unique_stop = find_xvals(*d.distributions[0].parameters)
        para_start, para_stop = find_xvals(*d.distributions[1].parameters)
        xvals = np.linspace(min(unique_start, para_start), max(unique_stop, para_stop), 1000)
        paravals = np.array([norm.pdf(x, *d.distributions[1].parameters) for x in xvals])
        uniquevals = np.array([norm.pdf(x, *d.distributions[0].parameters) for x in xvals])
        x = xvals[np.argwhere(np.diff(np.sign(uniquevals - paravals)) != 0).reshape(-1)[0]]
        return reverse_transform(x)

    biotype_df = pd.merge(aln_eval_df, ref_df, on=['TranscriptId'])
    r = []  # will hold the labels
    identity_cutoffs = []

    for biotype, df in biotype_df.groupby('TranscriptBiotype'):
        biotype_unique = df[df.Paralogy == 1]
        biotype_not_unique = df[df.Paralogy != 1]
        if len(biotype_not_unique) == 0:
            # No paralogous mappings implies all passing
            r.extend([[aln_id, 'Passing'] for aln_id in df.AlignmentId])
            identity_cutoffs.append([biotype, None])
        elif len(biotype_unique) == 0:
            # Only paralogous mappings implies all failing
            r.extend([[aln_id, 'Failing'] for aln_id in df.AlignmentId])
            identity_cutoffs.append([biotype, None])
        else:
            unique_transformed = transform_data(biotype_unique.TransMapIdentity)
            para_transformed = transform_data(biotype_not_unique.TransMapIdentity)
            combined = pd.concat([unique_transformed, para_transformed])
            unique_mu, unique_sigma = calculate_starting_values(unique_transformed)
            para_mu, para_sigma = calculate_starting_values(para_transformed)
            d = GeneralMixtureModel([NormalDistribution(unique_mu, unique_sigma),
                                     NormalDistribution(para_mu, para_sigma)])
            _ = d.fit(combined)
            cutoff = find_cutoff(d)
            identity_cutoffs.append([biotype, cutoff])
            labels = ['Passing' if ident >= cutoff else 'Failing' for ident in df.TransMapIdentity]
            r.extend([[aln_id, label] for aln_id, label in zip(*[df.AlignmentId, labels])])

    # turn r into a dataframe
    r_df = pd.DataFrame(r)
    r_df.columns = ['AlignmentId', 'TranscriptClass']
    ident_df = pd.DataFrame(identity_cutoffs)
    ident_df.columns = ['TranscriptBiotype', 'IdentityCutoff']
    return pd.merge(biotype_df, r_df, on='AlignmentId'), ident_df


def resolve_paralogs(updated_aln_eval_df):
    """
    Resolves paralogs based on likelihood to come from the two lognormal distributions.

    1. If only one paralog is more likely under the ortholog model, discard the others.
    2. If more than one paralog are more likely under the ortholog model, or we were unable to fit a model, resolve
       based on a heuristic combination of synteny score and alignment badness.
       score: 0.10 * coverage + 0.65 * (1 - badness) + 0.25 * (synteny / 6)

    :param updated_aln_eval_df: DataFrame produced by fit_distributions()
    :return: tuple of (metrics_dict, filtered DataFrame)
    """
    def score_aln(s):
        return 0.10 * s.TransMapCoverage + 0.65 * (1 - s.TransMapBadness) + 0.25 * (1.0 * s.Synteny / 6)

    def apply_label(s):
        return s.TranscriptClass if s.ParalogStatus != 'NotConfident' else 'Failing'

    updated_aln_eval_df['Score'] = updated_aln_eval_df.apply(score_aln, axis=1)
    updated_aln_eval_df = updated_aln_eval_df.sort_values(by='Score', ascending=False)

    paralog_status = []  # stores the results for a new column
    paralog_metrics = {biotype: {'Alignments discarded': 0, 'Transcripts resolved by model prediction': 0,
                                 'Transcripts resolved by synteny heuristic': 0, 'Transcripts arbitarily resolved': 0}
                       for biotype in set(updated_aln_eval_df.TranscriptBiotype)}

    for tx, df in updated_aln_eval_df.groupby('TranscriptId'):
        if len(df) == 1:  # no paralogs
            paralog_status.append([df.AlignmentId.iloc[0], None])
            continue
        biotype = df.TranscriptBiotype.iloc[0]
        passing = df[df.TranscriptClass == 'Passing']
        if len(passing) == 1:  # we can pick one passing member
            paralog_metrics[biotype]['Alignments discarded'] += 1
            paralog_metrics[biotype]['Transcripts resolved by model prediction'] += 1
            paralog_status.append([df.AlignmentId.iloc[0], 'Confident'])
        else:
            highest_score_df = df[df.Score == df.iloc[0].Score]
            if len(highest_score_df) == 1:
                paralog_metrics[biotype]['Transcripts resolved by synteny heuristic'] += 1
                paralog_metrics[biotype]['Alignments discarded'] += 1
                paralog_status.append([highest_score_df.AlignmentId.iloc[0], 'Confident'])
            else:
                paralog_metrics[biotype]['Transcripts arbitarily resolved'] += 1
                paralog_metrics[biotype]['Alignments discarded'] += 1
                paralog_status.append([highest_score_df.AlignmentId.iloc[0], 'NotConfident'])

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
    For each chromosome a gene maps to, find the weighted average of:
        0.15 * coverage + 0.65 * (1 - badness) + 0.25 * (synteny / 6)
    Pick the chromosome with the highest value of this metric.

    If the gene biotype is protein_coding, only count coding transcripts.

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
        """finds the sum of of (synteny + 1) * distance for each chrom in chroms"""
        r = {}
        tot = len(rec)
        for chrom, vals in chroms.iteritems():
            r[chrom] = (1.0 * len(vals) / tot) * np.mean(zip(*vals)[1])
        return r

    def score_aln(s):
        """scores an alignment based on badness, distance and synteny"""
        return 0.15 * s.TransMapCoverage + 0.65 * (1 - s.TransMapBadness) + 0.25 * (1.0 * s.Synteny / 6)

    def find_coding_chromosomes(rec):
        """create a mapping of genome sequence names to associated protein coding alignment ids and synteny scores"""
        chroms = collections.defaultdict(list)
        for _, s in rec.iterrows():
            if s.TranscriptBiotype != 'protein_coding':
                continue
            chroms[tx_dict[s.AlignmentId].chromosome].append([s.AlignmentId, score_aln(s)])
        return chroms

    def find_chromosomes(rec):
        """create a mapping of genome sequence names to associated alignment ids and synteny scores"""
        chroms = collections.defaultdict(list)
        for _, s in rec.iterrows():
            chroms[tx_dict[s.AlignmentId].chromosome].append([s.AlignmentId, score_aln(s)])
        return chroms

    def find_best_chroms(chroms):
        """finds the best chrom in the chroms dict based on the lowest chrom_metric score"""
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
            chroms = find_coding_chromosomes(rec)
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
        split_gene_metrics['Number of transcripts removed'] += 1
        other_contigs = ','.join(set(chroms.keys()) - {best_chroms[0]})
        split_status.extend([[tx_id, other_contigs] for tx_id in rec.TranscriptId])

    split_df = pd.DataFrame(split_status)
    split_df.columns = ['TranscriptId', 'GeneAlternateContigs']
    merged_df = pd.merge(paralog_filtered_df, split_df, on='TranscriptId')
    final_df = merged_df[~merged_df['AlignmentId'].isin(alignment_ids_to_remove)]
    return split_gene_metrics, final_df


def create_new_table(paralog_filtered_df):
    """
    Clean up this dataframe to write to SQL. I am not using the long form here because it is too hard to reload when
    the columns are not numeric.
    :param paralog_filtered_df: output from either resolve_split_genes() if flag set else resolve_paralogs()
    :return: dataframe to be written to sql
    """
    df = paralog_filtered_df[['GeneId', 'TranscriptId', 'AlignmentId',
                              'TranscriptClass', 'ParalogStatus', 'GeneAlternateContigs']]
    return df.set_index(['TranscriptId', 'AlignmentId'])

