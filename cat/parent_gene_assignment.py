"""
A set of functions to perform parental gene assignment in the AugustusPB/AugustusCGP modules
"""
import pandas as pd
import itertools
import collections
import tools.procOps
import tools.fileOps
import tools.mathOps
import tools.transcripts
import tools.intervals


def assign_parents(filtered_tm_gp, unfiltered_tm_gp, chrom_sizes, denovo_gp, min_distance=0.9, stranded=True):
    """
    Main function for assigning parental genes. Parental gene assignment methodology:
    A) Each denovo transcript is evaluated for overlapping any transMap transcripts.
    Overlap is defined as having at least 1 exonic base shared.
    B) If a denovo transcript is assigned to more than one gene, then this is attempted to be resolved.
    Resolution occurs by looking first at the transMap themselves. If any transMap projections overlap each other
    with a Jaccard metric > 0.001, then we call this as a badAnnotOrTm. These will be discarded unless all splices
    are supported.
    Next, we look at the asymmetric distance between this prediction and the gene intervals. If this difference
    in these distances is over min_distance for all comparisons, we call this rescued and it can be incorporated.
    Otherwise, this transcript is tagged ambiguousOrFusion.
    """
    filtered_transmap_dict = tools.transcripts.get_gene_pred_dict(filtered_tm_gp, stranded)
    unfiltered_transmap_dict = tools.transcripts.get_gene_pred_dict(unfiltered_tm_gp, stranded)
    filtered_ids = unfiltered_transmap_dict.keys() - filtered_transmap_dict.keys()

    tm_chrom_dict = create_chrom_dict(unfiltered_transmap_dict, chrom_sizes)
    denovo_dict = tools.transcripts.get_gene_pred_dict(denovo_gp, stranded)
    denovo_chrom_dict = create_chrom_dict(denovo_dict)

    # begin parent gene assignment
    r = []
    for chrom, tm_tx_by_chromosome in tm_chrom_dict.items():
        for denovo_tx_id, denovo_tx in denovo_chrom_dict[chrom].items():
            # find the names of both filtered and unfiltered transMap transcript IDs that overlap
            unfiltered_overlapping_tm_txs = find_tm_overlaps(denovo_tx, tm_tx_by_chromosome)
            filtered_overlapping_tm_txs = {tx for tx in unfiltered_overlapping_tm_txs if tx.name not in filtered_ids}
            # extract only gene names for the filtered set
            filtered_gene_ids = {tx.name2 for tx in filtered_overlapping_tm_txs}
            if len(filtered_gene_ids) > 1:  # we have more than one match, so resolve it
                resolved_name, resolution_method = resolve_multiple_genes(denovo_tx, filtered_overlapping_tm_txs,
                                                                          min_distance)
            elif len(filtered_gene_ids) == 1:  # yay, we have exactly one match
                resolved_name = list(filtered_gene_ids)[0]
                resolution_method = None
            else:
                resolved_name = resolution_method = None  # we have no matches, which means putative novel
            # find only genes for the unfiltered set that are not present in the filtered set
            alternative_gene_ids = {tx.name2 for tx in unfiltered_overlapping_tm_txs} - {resolved_name}
            alternative_gene_ids = ','.join(alternative_gene_ids) if len(alternative_gene_ids) > 0 else None
            r.append([denovo_tx.name, resolved_name, alternative_gene_ids, resolution_method])

    combined_alternatives = pd.DataFrame(r, columns=['TranscriptId', 'AssignedGeneId', 'AlternativeGeneIds',
                                                     'ResolutionMethod'])
    combined_alternatives = combined_alternatives.set_index('TranscriptId')
    return combined_alternatives


def create_chrom_dict(tx_dict, chrom_sizes=None):
    """
    Split up a dictionary of Transcript objects by chromosome. Add in extra chromosomes based on a sizes file
    """
    chrom_dict = collections.defaultdict(dict)
    for tx_id, tx in tx_dict.items():
        chrom_dict[tx.chromosome][tx_id] = tx
    if chrom_sizes is not None:
        for chrom, size in tools.fileOps.iter_lines(chrom_sizes):
            if chrom not in chrom_dict:
                chrom_dict[chrom] = {}
    return chrom_dict


def find_tm_overlaps(denovo_tx, tm_tx_dict):
    """Find overlap with transMap transcripts first on a genomic scale then an exonic scale"""
    r = []
    for tx in tm_tx_dict.values():
        if tx.interval.intersection(denovo_tx.interval) is not None:
            # make sure that we have exon overlap
            if ensure_exon_overlap(tx, denovo_tx) is True:
                r.append(tx)
    return r


def ensure_exon_overlap(tx, denovo_tx):
    """Do these two transcripts have at least 1 exonic base of overlap?"""
    for tm_exon in tx.exon_intervals:
        for denovo_exon in denovo_tx.exon_intervals:
            if tm_exon.overlap(denovo_exon) is True:
                return True
    return False


def resolve_multiple_genes(denovo_tx, overlapping_tm_txs, min_distance):
    """
    Resolve multiple assignments based on the following rules:
    """
    # use Jaccard metric to determine if the problem lies with transMap or annotation
    tm_txs_by_gene = tools.transcripts.group_transcripts_by_name2(overlapping_tm_txs)
    tm_jaccards = [find_highest_gene_jaccard(x, y) for x, y in itertools.combinations(list(tm_txs_by_gene.values()), 2)]
    if any(x > 0.001 for x in tm_jaccards):
        return None, 'badAnnotOrTm'
    # calculate asymmetric difference for this prediction
    scores = collections.defaultdict(list)
    for tx in overlapping_tm_txs:
        scores[tx.name2].append(calculate_asymmetric_closeness(denovo_tx, tx))
    best_scores = {gene_id: max(scores[gene_id]) for gene_id in scores}
    high_score = max(best_scores.values())
    if all(high_score - x >= min_distance for x in best_scores.values() if x != high_score):
        best = sorted(iter(best_scores.items()), key=lambda gene_id_score: gene_id_score[1])[-1][0]
        return best, 'rescued'
    else:
        return None, 'ambiguousOrFusion'


def find_highest_gene_jaccard(gene_list_a, gene_list_b):
    """
    Calculates the overall distance between two sets of transcripts by finding their distinct exonic intervals and then
    measuring the Jaccard distance.
    """
    def find_interval(gene_list):
        gene_intervals = set()
        for tx in gene_list:
            gene_intervals.update(tx.exon_intervals)
        gene_intervals = tools.intervals.gap_merge_intervals(gene_intervals, 0)
        return gene_intervals

    a_interval = find_interval(gene_list_a)
    b_interval = find_interval(gene_list_b)
    return tools.intervals.calculate_bed12_jaccard(a_interval, b_interval)


def calculate_asymmetric_closeness(denovo_tx, tm_tx):
    """
    Calculates the asymmetric closeness between two transcripts. This allows for denovo predictions that are subsets
    of existing annotations to get more weight.

    closeness = length(intersection) / length(denovo) in chromosome space
    """
    intersection = denovo_tx.interval.intersection(tm_tx.interval)
    if intersection is None:
        return 0
    return tools.mathOps.format_ratio(len(intersection), len(denovo_tx.interval))
