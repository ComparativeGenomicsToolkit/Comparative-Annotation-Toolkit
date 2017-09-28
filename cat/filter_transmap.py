"""
Filtering transMap.

globalNearBest is used to determine the single best alignment for a given transcript. If transcripts for a locus
end up disjoint, try to resolve this, rescuing lower-quality alignments as necessary.

"""
import json
import os
import logging
import itertools
import collections
import numpy as np
import pandas as pd
from networkx import Graph
from networkx.algorithms import connected_components
import tools.nameConversions
import tools.transcripts
import tools.psl
import tools.mathOps
import tools.procOps
import tools.fileOps
import tools.intervals
import tools.sqlInterface
from bx.intervals.cluster import ClusterTree

pd.options.mode.chained_assignment = None
logger = logging.getLogger(__name__)


def filter_transmap(tm_psl, ref_psl, tm_gp, genome, db_path, psl_tgt, minimum_paralog_coverage, local_near_best,
                    json_tgt):
    """
    Entry point for transMap filtering.
    :param tm_psl: input PSL
    :param ref_psl: reference fake PSL
    :param tm_gp: genePred from tm_psl
    :param genome: What genome are we working on? For logging purposes.
    :param db_path: Path to reference database, to get gene name to transcript name mapping
    :param psl_tgt: luigi.LocalTarget() object for PSL output
    :param minimum_paralog_coverage: Minimum coverage of a filtered alignment to be considered a paralogy
    :param local_near_best: localNearBest value to pass to PslCDnaFilter
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

    # find paralogies by using localNearBest
    with tools.fileOps.TemporaryFilePath() as local_tmp:
        cmd = [['sed', 's/\-[0-9]\+//', tmp_size_filtered],  # strip unique identifiers for comparative filters
               ['pslCDnaFilter', '-localNearBest={}'.format(local_near_best),
                '-minCover={}'.format(minimum_paralog_coverage / 100), '-verbose=0',
                '-minSpan=0.2', '/dev/stdin', '/dev/stdout']]
        tools.procOps.run_proc(cmd, stdout=local_tmp)
        paralogy_alns = list(tools.psl.psl_iterator(local_tmp))

    # load localBest IDs by using the hash table to figure out which ones we had
    local_best = {unfiltered[unfiltered_hash_table[hash_aln(aln)]] for aln in paralogy_alns}
    # report counts by biotype
    grouped = tools.psl.group_alignments_by_qname(local_best)
    metrics = {'Paralogy': collections.defaultdict(lambda: collections.Counter())}
    for tx_id, aln_ids in grouped.iteritems():
        biotype = transcript_biotype_map[tx_id]
        metrics['Paralogy'][biotype][len(aln_ids)] += 1

    # now perform globalNearBest to resolve exact orthologs
    with tools.fileOps.TemporaryFilePath() as stats_tmp, tools.fileOps.TemporaryFilePath() as filtered_psl:
        cmd = [['sed', 's/\-[0-9]\+//', tmp_size_filtered],  # strip unique identifiers for comparative filters
               ['pslCDnaFilter', '-globalNearBest=0', '-minCover=0.1', '-minSpan=0.2',
                '-statsOut={}'.format(stats_tmp), '/dev/stdin', '/dev/stdout']]
        tools.procOps.run_proc(cmd, stdout=filtered_psl)

        # load stats
        stats = parse_stats(stats_tmp)

        # load globalBest IDs by using the hash table to figure out which one we had
        filtered = tools.psl.get_alignment_dict(filtered_psl)
        best_ids = {unfiltered_hash_table[hash_aln(aln)] for aln in filtered.itervalues()}

    # run pslCDnaFilter again, with no options, to get scores
    with tools.fileOps.TemporaryFilePath() as tmp_verbose:
        cmd = ['pslCDnaFilter', '-verbose=5', tmp_size_filtered, '/dev/stdout']
        tools.procOps.run_proc(cmd, stderr=tmp_verbose, stdout='/dev/null')
        scores = parse_verbose(tmp_verbose)

    # resolve split genes using the scores and the best IDs
    resolved_df, num_rescued, split_gene_metrics = resolve_split_genes(best_ids, size_filtered,
                                                                       unfiltered_tx_dict, transcript_gene_map, scores)
    metrics['Split Genes'] = split_gene_metrics

    # add paralogy info
    resolved_df['Paralogy'] = [','.join({y.q_name for y in grouped[tx_id]} - {aln_id}) for tx_id, aln_id in
                               zip(*[resolved_df['TranscriptId'], resolved_df['AlignmentId']])]
    resolved_df['Paralogy'] = resolved_df.Paralogy.replace('', np.nan)

    # record stats and report to log
    stats['Rescued'] = num_rescued
    stats['Max Span Distance'] = num_too_long
    logger.info('Dropped {Coverage Filter:,} alignments due to low coverage, {Min Span Distance:,} '
                'alignments due to low spanning distance, '
                '{Max Span Distance:,} due to suspiciously large spanning distance, '
                'and {Paralog Filter:,} alignments during ortholog resolution. '
                'Split gene resolution rescued {Rescued:,} transcripts '
                'for a total of {Total:,} orthologs for genome {genome}.'.format(genome=genome, Total=len(resolved_df),
                                                                                 **stats))
    metrics['Orthology'] = stats

    # write the paralog resolved PSL
    with psl_tgt.open('w') as outf:
        for aln_id in resolved_df.AlignmentId:
            aln = unfiltered[aln_id]
            tools.fileOps.print_row(outf, aln.psl_string())

    # write the JSON
    tools.fileOps.ensure_file_dir(json_tgt.path)
    with json_tgt.open('w') as outf:
        json.dump(metrics, outf)

    os.remove(tmp_size_filtered)

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

    r = {}
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


def split_cluster(cluster):
    """If a cluster contains groups with no exonic overlap, split them up into new clusters"""
    # edge case -- if we have only 1 transcript, this algorithm fails
    if len(cluster) == 1:
        return [[0]]
    ct = ClusterTree(0, 1)
    for i, tx in enumerate(cluster):
        for e in tx.exon_intervals:
            ct.insert(e.start, e.stop, i)
    # now for each disjoint interval, which transcripts did we hit?
    indices = [set(indices) for start, end, indices in ct.getregions()]
    G = Graph()
    for s in indices:
        # this alignment is all alone
        if len(s) == 1:
            G.add_node(list(s)[0])
        for i, j in itertools.combinations(s, 2):
            G.add_edge(i, j)
    return list(connected_components(G))


def resolve_split_genes(best_ids, unfiltered, unfiltered_tx_dict, transcript_gene_map, scores):
    """
    Resolve split genes. Cluster all filtered transcripts together. If any genes are in different places,
    then pick the place with the highest sum of scores. Any transcripts in other locations get their other alignment
    with the highest score to be rescued.
    """
    # grab all best ID transcripts and cluster by-gene
    by_gene_best = collections.defaultdict(list)
    by_gene = collections.defaultdict(lambda: collections.defaultdict(list))
    for aln_id, aln in unfiltered.iteritems():
        tx_id = tools.nameConversions.strip_alignment_numbers(aln_id)
        gene_id = transcript_gene_map[tx_id]
        by_gene[gene_id][tx_id].append(aln)
        if aln_id in best_ids:
            by_gene_best[gene_id].append(aln)

    # keep track of genes which have to be resolved and if they are on the same contig or different contigs
    split_gene_metrics = {'Number of contig split genes': 0,
                          'Number of intra-contig split genes': 0}
    num_rescued = 0
    resolved = []
    # cluster the best alignments of each gene separately
    for gene_id, alns in by_gene_best.iteritems():

        # construct by-target cluster trees
        cts = collections.defaultdict(lambda: ClusterTree(0, 1))
        for i, aln in enumerate(alns):
            cts[aln.t_name].insert(aln.t_start, aln.t_end, i)

        if len(cts) > 1:  # we are split across contigs, record this for plots
            split_gene_metrics['Number of contig split genes'] += 1

        # now divide these cluster trees by ensuring exonic overlap
        scored_clusters = {}
        for chrom, ct in cts.iteritems():
            for start, end, aln_indices in ct.getregions():
                ct_alns = [alns[i] for i in aln_indices]
                ct_txs = [unfiltered_tx_dict[x.q_name] for x in ct_alns]
                for split_cluster_indices in split_cluster(ct_txs):
                    cluster_txs = [ct_txs[i] for i in split_cluster_indices]
                    cluster_alns = [ct_alns[i] for i in split_cluster_indices]
                    cluster_score = sum([scores[x.q_name] for x in cluster_alns])
                    sc_start = min(x.start for x in cluster_txs)
                    sc_end = max(x.stop for x in cluster_txs)
                    cluster_interval = tools.intervals.ChromosomeInterval(chrom, sc_start, sc_end, '.')
                    scored_clusters[cluster_interval] = [cluster_alns, cluster_score]

        if len(scored_clusters) == 1:
            # no split, just record and move on
            for aln in alns:
                tx_id = tools.nameConversions.strip_alignment_numbers(aln.q_name)
                resolved.append([gene_id, tx_id, aln.q_name, None])
            continue

        if len({x.chromosome for x in scored_clusters}) != len(scored_clusters):  # we are split in the same contig
            split_gene_metrics['Number of intra-contig split genes'] += 1

        # resolve by picking highest scoring clusters
        best_score = max([x[1] for x in scored_clusters.itervalues()])
        best_clusters = [[i, a] for i, (a, s) in scored_clusters.iteritems() if s == best_score]

        if len(best_clusters) > 1:
            # pick whichever cluster has the highest overall coverage
            i, best_cluster = sorted(best_clusters, key=lambda c: sum(x.coverage for x in c[1]))[-1]
        else:
            i, best_cluster = best_clusters[0]
        # keep map of tx_id to aln_id so we only end up with a single representative
        best_cluster_ids = {tools.nameConversions.strip_alignment_numbers(x.q_name): x.q_name for x in best_cluster}

        # construct a string marking the alternative loci
        alt_loci = scored_clusters.viewkeys() - {i}
        alt_loci = ','.join('{}:{}-{}'.format(x.chromosome, x.start, x.stop) for x in alt_loci)

        # construct the full set of exonic intervals in this cluster so that we don't grab alignments with no
        # exonic overlap. This hopefully prevents crazy transMaps from being pulled in.
        exons = [unfiltered_tx_dict[aln_id].exon_intervals for aln_id in best_cluster_ids.itervalues()]
        cluster_exons = tools.intervals.gap_merge_intervals(itertools.chain.from_iterable(exons), 0)

        # perform rescue process by finding all transcripts that overlap the given interval in the unfiltered set
        for tx_id, alns in by_gene[gene_id].iteritems():
            # find all alignments in the unfiltered set that overlap this locus
            rescue_candidates = []
            for aln in alns:
                aln_tx = unfiltered_tx_dict[aln.q_name]
                if any(not tools.intervals.interval_not_intersect_intervals(cluster_exons, e) for e in aln_tx.exon_intervals):
                    rescue_candidates.append(aln.q_name)
            if len(rescue_candidates) == 0:
                continue
            # use score to pick the best
            scored_candidates = [[aln_id, scores[aln_id]] for aln_id in rescue_candidates]
            best_candidate = sorted(scored_candidates, key=lambda x: x[1])[-1][0]
            best_cluster_ids[tools.nameConversions.strip_alignment_numbers(best_candidate)] = best_candidate
            num_rescued += 1

        for tx_id, aln_id in best_cluster_ids.iteritems():
            resolved.append([gene_id, tx_id, aln_id, alt_loci])
    resolved_df = pd.DataFrame(resolved, columns=['GeneId', 'TranscriptId', 'AlignmentId', 'GeneAlternateLoci'])
    return resolved_df, num_rescued, split_gene_metrics

