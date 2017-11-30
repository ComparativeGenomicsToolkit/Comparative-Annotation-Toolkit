"""
Filtering transMap.

"""
import os
import json
import logging
import collections
import pandas as pd
import tools.nameConversions
import tools.transcripts
import tools.psl
import tools.mathOps
import tools.procOps
import tools.fileOps
import tools.intervals
import tools.sqlInterface

pd.options.mode.chained_assignment = None
logger = logging.getLogger(__name__)


def filter_transmap(tm_psl, ref_psl, tm_gp, db_path, psl_tgt, global_near_best, json_tgt):
    """
    Entry point for transMap filtering.
    :param tm_psl: input PSL
    :param ref_psl: reference fake PSL
    :param tm_gp: genePred from tm_psl
    :param db_path: Path to reference database, to get gene name to transcript name mapping
    :param psl_tgt: luigi.LocalTarget() object for PSL output
    :param global_near_best: globalNearBest value to pass to PslCDnaFilter
    :param json_tgt: luigi.localTarget() object for JSON output
    :return:
    """
    # load all of the input alignments
    unfiltered = tools.psl.get_alignment_dict(tm_psl)
    unfiltered_tx_dict = tools.transcripts.get_gene_pred_dict(tm_gp)
    ref_psl_dict = tools.psl.get_alignment_dict(ref_psl)

    # pre-filter out suspiciously large spans
    size_filtered, num_too_long = ref_span(unfiltered, ref_psl_dict)
    tmp_size_filtered = tools.fileOps.get_tmp_file()
    with open(tmp_size_filtered, 'w') as outf:
        for aln in size_filtered.itervalues():
            tools.fileOps.print_row(outf, aln.psl_string())

    # get transcript -> gene map
    transcript_gene_map = tools.sqlInterface.get_transcript_gene_map(db_path)
    # get transcript -> biotype map for metrics
    transcript_biotype_map = tools.sqlInterface.get_transcript_biotype_map(db_path)
    # get annotation information for common names
    annotation_df = tools.sqlInterface.load_annotation(db_path)
    gene_name_map = dict(zip(annotation_df.GeneId, annotation_df.GeneName))

    # Construct a hash of alignment metrics to alignment IDs
    # The reason for this is that pslCDnaFilter rearranges them internally, so we lose order information

    def hash_aln(aln):
        """Hacky way to hash an alignment"""
        return hash(tuple([aln.t_name, aln.t_start, aln.t_end, aln.matches, aln.mismatches, aln.block_count,
                           tools.nameConversions.strip_alignment_numbers(aln.q_name),
                           tuple(aln.t_starts), tuple(aln.q_starts), tuple(aln.block_sizes)]))

    unfiltered_hash_table = {}
    for aln_id, aln in unfiltered.iteritems():
        unfiltered_hash_table[hash_aln(aln)] = aln_id
    assert len(unfiltered_hash_table) == len(unfiltered)

    with tools.fileOps.TemporaryFilePath() as local_tmp:
        cmd = [['sed', 's/\-[0-9]\+//', tmp_size_filtered],  # strip unique identifiers for comparative filters
               ['pslCDnaFilter', '-globalNearBest={}'.format(global_near_best),
                '-minCover=0.1', '-verbose=0',
                '-minSpan=0.2', '/dev/stdin', '/dev/stdout']]
        tools.procOps.run_proc(cmd, stdout=local_tmp)
        filtered_alns = list(tools.psl.psl_iterator(local_tmp))

    # load globalBest IDs by using the hash table to figure out which ones we had
    global_best = {unfiltered[unfiltered_hash_table[hash_aln(aln)]] for aln in filtered_alns}
    global_best_txs = [unfiltered_tx_dict[aln.q_name] for aln in global_best]

    # report counts by biotype
    grouped = tools.psl.group_alignments_by_qname(global_best)
    metrics = {'Paralogy': collections.defaultdict(lambda: collections.Counter())}
    paralogy_df = []
    for tx_id, alns in grouped.iteritems():
        biotype = transcript_biotype_map[tx_id]
        paralogy_df.append([tx_id, ','.join(sorted([x.q_name for x in alns]))])
        metrics['Paralogy'][biotype][len(alns)] += 1

    paralogy_df = pd.DataFrame(paralogy_df, columns=['TranscriptId', 'Paralogy'])

    # run pslCDnaFilter again, with no options, to get scores
    with tools.fileOps.TemporaryFilePath() as tmp_verbose:
        cmd = ['pslCDnaFilter', '-verbose=5', tmp_size_filtered, '/dev/stdout']
        tools.procOps.run_proc(cmd, stderr=tmp_verbose, stdout='/dev/null')
        scores = parse_verbose(tmp_verbose)

    # now coding and non-coding genes are split up. Coding genes are any genes who have a transcript with an ORF
    # the reason to do this is that running clusterGenes on full transcripts can lead to false fusions because
    # overlapping UTR intervals are real.

    # identify all genes that are non-coding. Non-coding is defined as genes who have no ORFs
    global_best_by_gene = tools.transcripts.group_transcripts_by_name2(global_best_txs)
    coding_genes = {gene_id for gene_id, tx_list in global_best_by_gene.iteritems()
                    if any(x.cds_size > 0 for x in tx_list)}

    with tools.fileOps.TemporaryFilePath() as coding_tmp, tools.fileOps.TemporaryFilePath() as noncoding_tmp, \
        tools.fileOps.TemporaryFilePath() as coding_clusters, tools.fileOps.TemporaryFilePath() as noncoding_clusters:
        with open(coding_clusters, 'w') as out_coding, open(noncoding_clusters, 'w') as out_noncoding:
            for tx in global_best_txs:
                if tx.name2 in coding_genes:
                    tools.fileOps.print_row(out_coding, tx.get_gene_pred())
                else:
                    tools.fileOps.print_row(out_noncoding, tx.get_gene_pred())
        cmd = ['clusterGenes', '-cds', coding_tmp, 'no', coding_clusters]
        tools.procOps.run_proc(cmd)
        cmd = ['clusterGenes', noncoding_tmp, 'no', noncoding_clusters]
        tools.procOps.run_proc(cmd)
        coding_clustered = pd.read_csv(coding_tmp, sep='\t')
        noncoding_clustered = pd.read_csv(noncoding_tmp, sep='\t')

    coding_merged_df, coding_collapse_filtered = filter_clusters(coding_clustered, transcript_gene_map,
                                                                 gene_name_map, scores)
    noncoding_merged_df, noncoding_collapse_filtered = filter_clusters(noncoding_clustered, transcript_gene_map,
                                                                       gene_name_map, scores)

    merged_collapse_filtered = pd.concat([coding_collapse_filtered, noncoding_collapse_filtered])
    merged_df = pd.concat([coding_merged_df, noncoding_merged_df])

    # Now that these have been processed separately, two things must happen:
    # 1) All non-coding isoforms of coding genes must be re-added
    # 2) Alignments filtered by globalNearBest for a filtered gene should be rescued if they have sufficient coverage

    # first, group the putative rescue transcripts by gene ID. Require that they exist in scores because otherwise
    # that means they are weirdly overlapping
    high_cov_ids = {x.q_name for x in unfiltered.itervalues() if x.coverage > 0.5 and x.q_name in scores}
    high_cov_ids -= set(merged_collapse_filtered.gene)  # gene is alignment ID
    putative_rescue_txs = {tx for aln_id, tx in unfiltered_tx_dict.iteritems() if aln_id in high_cov_ids}
    unfiltered_by_gene = tools.transcripts.group_transcripts_by_name2(putative_rescue_txs)

    rescued_txs = []
    # for each gene ID that survived filtering, find their interval
    for gene_id, group in merged_collapse_filtered.groupby('gene_id'):
        assert len(set(group['#cluster'])) == 1
        tx_intervals = []
        for _, s in group.iterrows():
            tx_intervals.append(tools.intervals.ChromosomeInterval(s.chrom, s.txStart, s.txEnd, s.strand))
        tx_intervals = tools.intervals.hull_of_intervals(tx_intervals)
        assert tx_intervals is not None
        gene_interval = tx_intervals[0]
        for tx in unfiltered_by_gene[gene_id]:
            if tx.interval.overlap(gene_interval):
                rescued_txs.append(tx.name)

    # the final step is filtering for duplicates. Duplicates here means that we have multiple transMap
    # mapping to the same locus. Pick the highest scores
    combined_txs = rescued_txs + list(merged_collapse_filtered.gene)
    combined_tx_df = pd.DataFrame(combined_txs, columns=['AlignmentId'])
    combined_tx_df['score'] = [scores[x] for x in combined_tx_df.AlignmentId]
    combined_tx_df['TranscriptId'] = [tools.nameConversions.strip_alignment_numbers(x) for x in combined_tx_df.AlignmentId]
    combined_tx_df['GeneId'] = [transcript_gene_map[x] for x in combined_tx_df.TranscriptId]
    combined_tx_df = combined_tx_df.sort_values('score')
    combined_tx_df = combined_tx_df.groupby('TranscriptId', as_index=False).first()

    # construct the output DataFrame
    resolved_df = combined_tx_df.merge(merged_df, on='GeneId', how='left')
    resolved_df = resolved_df.drop('score', axis=1)

    # write the paralog resolved PSL
    with psl_tgt.open('w') as outf:
        for aln_id in resolved_df.AlignmentId:
            aln = unfiltered[aln_id]
            tools.fileOps.print_row(outf, aln.psl_string())

    # resolve split genes using the scores and the best IDs
    resolved_df, split_gene_metrics = resolve_split_genes(tmp_size_filtered, transcript_gene_map,
                                                          resolved_df, unfiltered_tx_dict)
    # add in paralogy calls from before
    resolved_df = resolved_df.merge(paralogy_df, on='TranscriptId')
    metrics['Split Genes'] = split_gene_metrics

    os.remove(tmp_size_filtered)

    # write the JSON
    tools.fileOps.ensure_file_dir(json_tgt.path)
    with json_tgt.open('w') as outf:
        json.dump(metrics, outf)

    return resolved_df.set_index(['GeneId', 'TranscriptId'])


def ref_span(aln_dict, ref_aln_dict, max_span=5):
    """
    --minSpan in pslCDnaFilter can lead to problematic results when transMap produces ultra-long transcripts.

    We reduce this problem by introducing a ref_span feature that looks at the genomic size of the reference transcripts
    and only keep those with a span up to max_span more than the reference genome.
    """
    # group by name
    grouped = collections.defaultdict(list)
    for aln_id, aln in aln_dict.iteritems():
        grouped[tools.nameConversions.strip_alignment_numbers(aln_id)].append(aln)

    r = collections.OrderedDict()
    for tx_id, aln_list in grouped.iteritems():
        ref_aln = ref_aln_dict[tx_id]
        ref_size = ref_aln.t_end - ref_aln.t_start
        ref_cutoff = ref_size * max_span
        alns = [aln for aln in aln_list if aln.t_end - aln.t_start <= ref_cutoff]
        for aln in alns:
            r[aln.q_name] = aln
    return r, len(aln_dict) - len(r)


def parse_stats(stats):
    """Parse the stats output, provide summary statistics to log"""
    stats = pd.read_csv(stats, sep='\t', names=['mode', 'seqs', 'alns'], index_col=0)
    # munge the stats and report them
    stats.index = [x.replace(' ', '') for x in stats.index]
    stats = stats.T
    stats_dict = {}
    if 'dropminCover:' in stats:
        stats_dict['Coverage Filter'] = int(stats['dropminCover:'].alns)
    else:
        stats_dict['Coverage Filter'] = 0
    if 'dropminSpan:' in stats:
        stats_dict['Min Span Distance'] = int(stats['dropminSpan:'].alns)
    else:
        stats_dict['Min Span Distance'] = 0
    if 'dropglobalBest:' in stats:
        stats_dict['Paralog Filter'] = int(stats['dropglobalBest:'].alns)
    else:
        stats_dict['Paralog Filter'] = 0
    return stats_dict


def parse_verbose(verbose):
    """Parse the verbose output to retain score information for resolution"""
    scores = {}
    for l in open(verbose):
        if l.startswith('align'):
            l = l.split()
            aln_id = l[-3].split(':')[0].split(']')[1]
            score = l[5]
            score = float(score.split('=')[1])
            scores[aln_id] = score
    return scores


def find_best_group(group, key):
    """
    Resolve a cluster by finding the highest average score. Key determines if we are currently resolving
    cluster or gene_id
    """
    avg_scores = group[[key, 'scores']].groupby(key, as_index=False).mean()
    return avg_scores.sort_values('scores', ascending=False).iloc[0][key]


def construct_alt_loci(group, best_cluster):
    """
    For paralogous genes, find the locations of alt loci
    """
    intervals = collections.defaultdict(list)
    for cluster_id, x in group.set_index('#cluster').iterrows():
        if cluster_id != best_cluster:
            intervals[x.chrom].append(tools.intervals.ChromosomeInterval(x.chrom, x.txStart, x.txEnd, '.'))
    merged_intervals = []
    for chrom, i in intervals.iteritems():
        merged_intervals.extend(tools.intervals.gap_merge_intervals(i, 1000))
    return ','.join('{}:{}-{}'.format(x.chromosome, x.start, x.stop) for x in merged_intervals)


def filter_clusters(clustered, transcript_gene_map, gene_name_map, scores):
    """
    Wrapper for taking the output of clusterGenes and filtering it
    """
    # add gene IDs and scores. clustered.gene is actually AlignmentId fields
    clustered['gene_id'] = [transcript_gene_map[tools.nameConversions.strip_alignment_numbers(x)] for x in clustered.gene]
    clustered['scores'] = [scores[x] for x in clustered.gene]

    to_remove = set()  # set of cluster IDs to remove
    alt_loci = []  # will become a DataFrame of alt loci to populate that field
    # any gene IDs with multiple clusters need to be resolved to resolve paralogies
    for gene_id, group in clustered.groupby('gene_id'):
        if len(set(group['#cluster'])) > 1:
            # pick the highest average scoring cluster
            best_cluster = find_best_group(group, '#cluster')
            best_cluster = int(best_cluster)
            alt_loci.append([gene_id, construct_alt_loci(group, best_cluster)])
            to_remove.update(set(group['#cluster']) - {best_cluster})
    paralog_filtered = clustered[~clustered['#cluster'].isin(to_remove)]
    paralog_df = pd.DataFrame(alt_loci, columns=['GeneId', 'GeneAlternateLoci'])

    # group by cluster ID to identify gene family collapse
    genes_to_remove = set()  # set of gene IDs to collapse
    collapsed_genes = []  # will become a DataFrame of collapsed genes
    for cluster_id, group in paralog_filtered.groupby('#cluster'):
        if len(set(group['gene_id'])) > 1:
            best_gene = find_best_group(group, 'gene_id')
            collapsed_gene_ids = set(group.gene_id) - {best_gene}
            collapsed_gene_names = {gene_name_map[x] for x in collapsed_gene_ids}
            to_remove.update(collapsed_gene_ids)
            collapsed_genes.append([best_gene, ','.join(collapsed_gene_ids), ','.join(collapsed_gene_names)])
    collapse_filtered = paralog_filtered[~paralog_filtered['gene_id'].isin(genes_to_remove)]
    collapsed_df = pd.DataFrame(collapsed_genes, columns=['GeneId', 'CollapsedGeneIds', 'CollapsedGeneNames'])
    merged_df = collapsed_df.merge(paralog_df, how='outer', on='GeneId')
    return merged_df, collapse_filtered


def find_split_genes(gene_id, g, resolved_interval, split_gene_data):
    """
    Determines if a group of alignments filtered by localNearBest contains putative split genes
    """
    intervals = collections.defaultdict(list)
    for aln in g:
        ref_i = tools.intervals.ChromosomeInterval(tools.nameConversions.strip_alignment_numbers(aln.q_name),
                                                   aln.q_start, aln.q_end, '.')
        tgt_i = tools.intervals.ChromosomeInterval(aln.t_name, aln.t_start, aln.t_end, aln.strand)
        intervals[ref_i].append(tgt_i)
    merged_intervals = tools.intervals.union_of_intervals(intervals.keys())
    if len(merged_intervals) > 1:
        alt_intervals = collections.defaultdict(list)
        for interval_list in intervals.itervalues():
            for i in interval_list:
                if not i.overlap(resolved_interval):
                    alt_intervals[i.chromosome].append(i)
        # merge by chromosome
        r = []
        for chrom, interval_list in alt_intervals.iteritems():
            r.extend(tools.intervals.gap_merge_intervals(interval_list, 0))
        # write metrics
        if len(alt_intervals) == 1 and alt_intervals.keys()[0] == resolved_interval.chromosome:
            split_gene_data['intra'].add(gene_id)
        else:
            split_gene_data['contig'].add(gene_id)
        if len(r) == 0:
            return None
        return ','.join(['{}:{}-{}'.format(i.chromosome, i.start, i.stop) for i in r])
    else:
        return None


def resolve_split_genes(tmp_size_filtered, transcript_gene_map, resolved_df, unfiltered_tx_dict):
    """
    Use localNearBest algorithm to determine split genes and populate that field
    """
    with tools.fileOps.TemporaryFilePath() as local_tmp:
        cmd = [['sed', 's/\-[0-9]\+//', tmp_size_filtered],  # strip unique identifiers for comparative filters
               ['pslCDnaFilter', '-localNearBest=0.05',
                '-minCover=0.1', '-verbose=0',
                '-minSpan=0.2', '/dev/stdin', '/dev/stdout']]
        tools.procOps.run_proc(cmd, stdout=local_tmp)
        filtered_alns = list(tools.psl.psl_iterator(local_tmp))

    # remove alignments that we didn't resolve
    resolved_ids = set(resolved_df.TranscriptId)
    filtered_alns = [x for x in filtered_alns if x.q_name in resolved_ids]
    grouped = tools.psl.group_alignments_by_qname(filtered_alns)

    # construct the transcript interval for resolved transcripts
    tx_intervals = {tx_id: unfiltered_tx_dict[aln_id].interval for
                    tx_id, aln_id in zip(resolved_df.TranscriptId, resolved_df.AlignmentId)}

    split_r = []
    # keep track of transcripts which have to be resolved and if they are on the same contig or different contigs
    split_gene_data = {'contig': set(), 'intra': set()}
    for tx_id, g in grouped.iteritems():
        gene_id = transcript_gene_map[tx_id]
        split_r.append([tx_id, find_split_genes(gene_id, g, tx_intervals[tx_id], split_gene_data)])
    split_df = pd.DataFrame(split_r, columns=['TranscriptId', 'PossibleSplitGeneLocations'])
    merged = split_df.merge(resolved_df, on='TranscriptId')

    # calculate the number of genes for metrics
    split_gene_metrics = {'Number of contig split genes': len(split_gene_data['contig']),
                          'Number of intra-contig split genes': len(split_gene_data['intra'])}

    return merged, split_gene_metrics
