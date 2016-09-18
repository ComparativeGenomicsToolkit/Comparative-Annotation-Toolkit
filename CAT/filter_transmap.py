"""
Filter transMap transcripts based on classification. These filtered transcripts are used as input to Augustus pipelines
as well as for the consensus gene set.

Filters:
1. LongTranscript. If any of the transcripts failed this classifier, remove it.
2. ResolveSplitGenes. If the flag to resolve split genes is set, resolve them based on synteny scores.
3. Paralogy. Attempt to resolve paralogs. If they are not resolvable, keep them. They will be passed to Augustus as
one transcript.
"""
import collections
import tools.sqlInterface
import tools.transcripts
import tools.nameConversions
import pandas as pd


def filter_transmap(filter_tm_args, out_target):
    """
    Entry point for transMap filtering.
    :param filter_tm_args: argparse Namespace produced by FilterTransMap.get_args()
    :param out_target: luigi.LocalTarget where the results will be written
    """
    # load database tables
    ref_df = tools.sqlInterface.load_reference(filter_tm_args.ref_db_path)
    aln_eval_df = tools.sqlInterface.load_alignment_evaluation(filter_tm_args.db_path)
    tx_dict = tools.transcripts.get_gene_pred_dict(filter_tm_args.tm_gp)

    # store metrics on filtering for plots. This will be written to disk for final plots
    metrics = {}

    # remove long transcripts
    long_ids_to_remove = filter_long_transcripts(aln_eval_df)
    metrics['Long Transcripts'] = len(long_ids_to_remove)
    aln_eval_df = filter_df(aln_eval_df, long_ids_to_remove)

    # remove super low coverage
    low_cov_ids_to_remove = filter_coverage(aln_eval_df)
    metrics['Coverage Filter'] = len(low_cov_ids_to_remove)
    aln_eval_df = filter_df(aln_eval_df, low_cov_ids_to_remove)

    # resolve paralogs
    paralog_metrics, paralog_ids_to_remove = resolve_paralogs(aln_eval_df)
    metrics['Paralogy'] = paralog_metrics
    aln_eval_df = filter_df(aln_eval_df, paralog_ids_to_remove)

    # resolve split genes, if user requested
    if filter_tm_args.resolve_split_genes is True:
        split_gene_metrics, split_ids_to_remove = resolve_split_genes(aln_eval_df, ref_df, tx_dict)
        metrics['Split Genes'] = split_gene_metrics
        aln_eval_df = filter_df(aln_eval_df, split_ids_to_remove)

    # write out the filtered transMap results
    with out_target.open('w') as outf:
        for aln_id in aln_eval_df.AlignmentId:
            tx = tx_dict[aln_id]
            outf.write('\t'.join(tx.get_gene_pred()) + '\n')


def filter_df(aln_eval_df, ids_to_remove):
    """remove IDs from aln_eval_df between each step"""
    return aln_eval_df[~aln_eval_df['AlignmentId'].isin(ids_to_remove)]


def filter_long_transcripts(aln_eval_df):
    """apply long transcript filter"""
    long_tx_df = aln_eval_df[aln_eval_df['LongTranscript'] == 1]
    return set(long_tx_df.AlignmentId)


def filter_coverage(aln_eval_df, cov_cutoff=0.1):
    """apply a very weak coverage filter. This reduces the number of very unlikely alignments"""
    low_cov_df = aln_eval_df[aln_eval_df['TransMapCoverage'] <= cov_cutoff]
    return set(low_cov_df.AlignmentId)


def resolve_paralogs(aln_eval_df):
    """
    Resolve paralogs based on synteny scores.
    1: If only one paralog has a non-zero synteny score, discard the others (easy case of retroposed copy)
    2: If one paralog for a given source transcript has a synteny score 3 or more than any other paralog, discard the
    others. (one region is a much better fit than the others)
    If 1/2 did not resolve this paralog, then we need to get to tricky ordering-based resolution. This will be done
    by AugustusTM(R) within that pipeline.

    :param aln_eval_df: DataFrame produced by load_alignment_evaluation
    :return: tuple of (metrics_dict, set of alignment ids to remove)
    """
    def pick_top(df, category):
        """throw away all alignment IDs not at the top of the sorted DataFrame"""
        alignment_ids = set(df.AlignmentId[1:])
        alignment_ids_to_remove.update(alignment_ids)
        paralog_metrics['Alignments discarded'] += len(alignment_ids)
        paralog_metrics[category] += 1

    alignment_ids_to_remove = set()
    paralog_metrics = collections.Counter()
    for tx, df in aln_eval_df.groupby('TranscriptId'):
        if len(df) == 1:  # no paralogs
            continue
        df = df.sort_values(by='Synteny', ascending=False)
        if df.iloc[0].Synteny - 3 > df.iloc[1].Synteny:
            pick_top(df, 'Transcripts resolved due to synteny delta')
        elif sum(df.Synteny) == df.iloc[0].Synteny:
            pick_top(df, 'Transcripts resolved due to one non-zero')
        else:
            paralog_metrics['Transcripts not resolved'] += 1
    return paralog_metrics, alignment_ids_to_remove


def resolve_split_genes(aln_eval_df, ref_df, tx_dict):
    """
    Resolves cases where transMap mapped a gene to different chromosomes. This is a useful feature to turn on
    if you have a high quality assembly, but may be problematic for highly fragmented assemblies.

    For each gene, if transcripts on that gene are on multiple sequences, a consensus finding process is performed.
    For each transMap transcript associated with this gene, sum up the synteny scores. Whichever sequence has the
    highest score wins. If there is a tie (or both are 0), pick the sequence with the lowest sum of badness.

    If the gene biotype is protein_coding, only count coding transcripts.

    :param aln_eval_df: DataFrame produced by load_alignment_evaluation
    :param ref_df: DataFrame produced by sqlInterface.load_reference()
    :param tx_dict: Dictionary mapping alignment IDs to GenePredTranscript objects for all alignment modes. We will
                    ignore the non-transMap ones. This is just to get the chromosomes.
    :return: tuple of (metrics_dict, set of alignment ids to remove)
    """
    def extract_gene_biotype(rec):
        """there should only be one gene biotype"""
        gene_biotypes = list(set(rec.GeneBiotype))
        assert len(gene_biotypes) == 1
        return gene_biotypes[0]

    def find_coding_chromosomes(rec):
        """create a mapping of genome sequence names to associated protein coding alignment ids and synteny scores"""
        chroms = collections.defaultdict(list)
        for aln_id, tx_biotype, synteny_score in zip(*[rec.AlignmentId, rec.TranscriptBiotype, rec.Synteny]):
            if tx_biotype != 'protein_coding':
                continue
            chroms[tx_dict[aln_id].chromosome].append([aln_id, synteny_score])
        return chroms

    def find_chromosomes(rec):
        """create a mapping of genome sequence names to associated alignment ids and synteny scores"""
        chroms = collections.defaultdict(list)
        for aln_id, synteny_score in zip(*[rec.AlignmentId, rec.Synteny]):
            chroms[tx_dict[aln_id].chromosome].append([aln_id, synteny_score])
        return chroms

    def synteny_metric(chroms):
        """sum up the synteny scores for each chromosome"""
        synteny_scores = collections.Counter()
        for chrom in chroms:
            for aln_id, synteny_score in chroms[chrom]:
                synteny_scores[chrom] += synteny_score
        return synteny_scores

    def find_best_chrom(synteny_scores):
        """find the best chromosome(s) based on synteny"""
        max_score = max(synteny_scores.values())
        best_chroms = [chrom for chrom, score in synteny_scores.iteritems() if score == max_score]
        return best_chroms

    def find_lowest_badness(chroms, best_chroms, rec):
        """we have a tie, pick based on the lowest sum of badness"""
        # filter transcript lists for best chrom only
        chrom_map = {chrom: {aln_id for aln_id, syn_score in chroms[chrom]} for chrom in best_chroms}
        # index the record for access
        rec = rec.set_index('AlignmentId')
        # sum up the badness scores
        badness_map = {}
        for chrom, aln_ids in chrom_map.iteritems():
            badness_map[chrom] = round(sum([rec.ix[aln_id].TransMapBadness for aln_id in aln_ids]), 5)
        min_badness = min([x for x in badness_map.itervalues()])
        best_chroms = [chrom for chrom, badness in badness_map.iteritems() if badness == min_badness]
        return best_chroms

    def find_names_to_remove(chrom_or_chroms, chroms):
        """return the list of names to remove as well as the count based on chrom_or_chroms"""
        if isinstance(chrom_or_chroms, str):
            chrom_list = [chrom_or_chroms]
        else:
            chrom_list = chrom_or_chroms
        aln_ids = set()
        for chrom in chrom_list:
            aln_ids.update([aln_id for aln_id, synteny_score in chroms[chrom]])
        return len(aln_ids), aln_ids

    # merge the aln_eval_df and ref_df, producing a gene-level name mapping
    merged_df = pd.merge(aln_eval_df, ref_df, on='TranscriptId', how='inner')

    # begin iteration, keeping track of metrics
    ids_to_remove = set()
    split_gene_metrics = collections.Counter()
    for gene, rec in merged_df.groupby('GeneId'):
        gene_biotype = extract_gene_biotype(rec)
        if gene_biotype == 'protein_coding':
            chroms = find_coding_chromosomes(rec)
            if len(chroms) == 0:  # there may be an edge case where a protein coding gene has no coding transcripts
                chroms = find_chromosomes(rec)
        else:
            chroms = find_chromosomes(rec)
        if len(chroms) == 1:  # no ambiguous chromosomes here
            continue
        split_gene_metrics['Number of genes affected'] += 1
        # try to resolve by synteny
        synteny_scores = synteny_metric(chroms)
        best_chroms = find_best_chrom(synteny_scores)
        if len(best_chroms) == 1:  # we have a winner, remove the rest
            split_gene_metrics['Resolved by synteny'] += 1
            n, i = find_names_to_remove(best_chroms[0], chroms)
            ids_to_remove.update(i)
        else:  # try to resolve by lowest sum of badness only on tie chromosomes
            best_chroms = find_lowest_badness(chroms, best_chroms, rec)
            if len(best_chroms) == 1:
                split_gene_metrics['Resolved by badness'] += 1
                n, i = find_names_to_remove(best_chroms[0], chroms)
                ids_to_remove.update(i)
            else:  # all hope is lost
                split_gene_metrics['Genes discarded'] += 1
                n, i = find_names_to_remove(chroms.keys(), chroms)
                ids_to_remove.update(i)
    return split_gene_metrics, ids_to_remove
