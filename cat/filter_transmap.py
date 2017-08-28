"""
Filtering transMap.

globalNearBest is used to determine the single best alignment for a given transcript. If transcripts for a locus
end up disjoint, try to resolve this, rescuing lower-quality alignments as necessary.

"""
import json
import logging
import collections
import numpy as np
import pandas as pd
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


def filter_transmap(tm_psl, genome, db_path, psl_tgt, minimum_paralog_coverage, local_near_best, json_tgt):
    """
    Entry point for transMap filtering.
    :param tm_psl: input PSL
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
        cmd = [['sed', 's/\-[0-9]\+//', tm_psl],  # strip unique identifiers for comparative filters
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
        cmd = [['sed', 's/\-[0-9]\+//', tm_psl],  # strip unique identifiers for comparative filters
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
        cmd = ['pslCDnaFilter', '-verbose=5', tm_psl, '/dev/stdout']
        tools.procOps.run_proc(cmd, stderr=tmp_verbose, stdout='/dev/null')
        scores = parse_verbose(tmp_verbose)

    resolved_df, to_remove, split_gene_metrics = resolve_split_genes(best_ids, unfiltered, transcript_gene_map, scores)
    metrics['Split Genes'] = split_gene_metrics

    # flatten these out to a dataframe for SQL
    r = []
    for tx_id, aln in filtered.iteritems():
        aln_id = unfiltered_hash_table[hash_aln(aln)]
        gene_id = transcript_gene_map[tx_id]
        r.append([gene_id, tx_id, aln_id])
    df = pd.DataFrame(r, columns=['GeneId', 'TranscriptId', 'AlignmentId'])
    df = df[~df.TranscriptId.isin(to_remove)]
    m = pd.concat([df, resolved_df])
    # add paralogy info
    m['Paralogy'] = [','.join({y.q_name for y in grouped[tx_id]} - {aln_id}) for tx_id, aln_id in
                     zip(*[m['TranscriptId'], m['AlignmentId']])]
    m['Paralogy'] = m.Paralogy.replace('', np.nan)

    # record stats and report to log
    stats['Rescued'] = len(m) - len(filtered) + len(to_remove)
    logger.info('Dropped {Coverage Filter:,} alignments due to low coverage, {Span Distance:,} '
                'alignments due to low spanning distance '
                'and {Paralog Filter:,} alignments during ortholog resolution. '
                'Split gene resolution rescued {Rescued:,} transcripts '
                'for a total of {Total:,} orthologs for genome {genome}.'.format(genome=genome, Total=len(m), **stats))
    metrics['Orthology'] = stats

    # write the paralog resolved PSL
    with psl_tgt.open('w') as outf:
        for aln_id in m.AlignmentId:
            aln = unfiltered[aln_id]
            tools.fileOps.print_row(outf, aln.psl_string())

    # write the JSON
    tools.fileOps.ensure_file_dir(json_tgt.path)
    with json_tgt.open('w') as outf:
        json.dump(metrics, outf)

    return m.set_index(['GeneId', 'TranscriptId'])


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
        stats_dict['Span Distance'] = int(stats['dropminSpan:'].alns)
    else:
        stats_dict['Span Distance'] = 0
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


def resolve_split_genes(best_ids, unfiltered, transcript_gene_map, scores):
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
    resolved = []
    # set of alignments in the globalNearBest set we want to remove after resolution
    to_remove = set()
    # cluster the best alignments of each gene separately
    for gene_id, alns in by_gene_best.iteritems():
        # construct by-target cluster trees
        cts = collections.defaultdict(lambda: ClusterTree(0, 1))
        for i, aln in enumerate(alns):
            cts[aln.t_name].insert(aln.t_start, aln.t_end, i)
        if len(cts) > 1:  # we are split across contigs, record this for plots
            split_gene_metrics['Number of contig split genes'] += 1
        # find the highest scoring cluster
        scored_clusters = {}
        for chrom, ct in cts.iteritems():
            for start, end, aln_idx in ct.getregions():
                interval = tools.intervals.ChromosomeInterval(chrom, start, end, '.')
                ct_alns = [alns[i] for i in aln_idx]
                cluster_score = np.mean([scores[x.q_name] for x in ct_alns])
                scored_clusters[interval] = [ct_alns, cluster_score]
        if len(scored_clusters) > 1:
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
            # rescue transcripts
            rescued_alns = []
            for tx_id, unfiltered_alns in by_gene[gene_id].iteritems():
                overlapping_alns = []
                for aln in unfiltered_alns:
                    if aln not in best_cluster:
                        aln_i = tools.intervals.ChromosomeInterval(aln.t_name, aln.t_start, aln.t_end, '.')
                        if aln_i.overlap(i):
                            overlapping_alns.append([aln, scores[aln.q_name]])
                if len(overlapping_alns) > 0:
                    best_overlap = sorted(overlapping_alns, key=lambda x: x[1])[-1][0]
                    rescued_alns.append(best_overlap)
                alt_contigs = ','.join({a.t_name for a in rescued_alns} - {aln.t_name})
                resolved.append([gene_id, tx_id, aln.q_name, alt_contigs])
                to_remove.add(tx_id)  # remove these tx IDs from globalNearBest set
    resolved_df = pd.DataFrame(resolved, columns=['GeneId', 'TranscriptId', 'AlignmentId', 'GeneAlternateContigs'])
    return resolved_df, to_remove, split_gene_metrics

