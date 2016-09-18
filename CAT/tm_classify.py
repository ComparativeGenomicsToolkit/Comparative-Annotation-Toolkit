"""
Classify transMap transcripts producing the TransMapEvaluation table for each genome's database

1. Paralogy: The # of times this transcript was aligned
2. AlnExtendsOffConfig: Does this alignment run off the end of a contig?
3. AlignmentPartialMap: Did this transcript not map completely?
4. AlnAbutsUnknownBases: Does this alignment have Ns immediately touching any exons?
5. AlnContainsUnknownBases: Are there any Ns within the transcript alignment?
6. LongAlignment: Did this transcript align in a insanely long fashion? Indicative of paralogy problems.
7. Synteny: If this transcript aligned more than once, assign a boolean based on synteny to whether this is the
    most probable transcript. This is used to filter for pseudogenes.
"""
import bisect
import collections
import pandas as pd
import tools.psl
import tools.transcripts
import tools.nameConversions
import tools.bio


# hard coded variables
# hard coded long transMap size. Bigger than 3 megabases is probably a spurious alignment.
long_tx_size = 3 * 10 ** 6


def tm_classify(tm_eval_args):
    """
    Runs alignment classification based on transMap PSLs, genePreds and the genome FASTA.
    :param tm_eval_args: argparse Namespace produced by EvaluateTransMap.get_args()
    :return: DataFrame
    """
    psl_dict = tools.psl.get_alignment_dict(tm_eval_args.tm_psl)
    gp_dict = tools.transcripts.get_gene_pred_dict(tm_eval_args.tm_gp)
    ref_gp_dict = tools.transcripts.get_gene_pred_dict(tm_eval_args.annotation_gp)
    fasta = tools.bio.get_sequence_dict(tm_eval_args.fasta)
    r = []
    paralog_count = paralogy(psl_dict)  # we have to count paralogs globally
    synteny_scores = synteny(ref_gp_dict, gp_dict)  # we also have to score synteny globally
    for aln_id, tx in gp_dict.iteritems():
        aln = psl_dict[aln_id]
        tx_id = tools.nameConversions.strip_alignment_numbers(aln_id)
        r.append([aln_id, tx_id, 'Paralogy', paralog_count[aln_id]])
        r.append([aln_id, tx_id, 'Synteny', synteny_scores[aln_id]])
        r.append([aln_id, tx_id, 'AlnExtendsOffContig', aln_extends_off_contig(aln)])
        r.append([aln_id, tx_id, 'AlnPartialMap', alignment_partial_map(aln)])
        r.append([aln_id, tx_id, 'AlnAbutsUnknownBases', aln_abuts_unknown_bases(tx, fasta)])
        r.append([aln_id, tx_id, 'AlnContainsUnknownBases', aln_contains_unknown_bases(tx, fasta)])
        r.append([aln_id, tx_id, 'LongTranscript', long_transcript(tx)])
        r.append([aln_id, tx_id, 'TransMapCoverage', aln.coverage])
        r.append([aln_id, tx_id, 'TransMapIdentity', aln.identity])
        r.append([aln_id, tx_id, 'TransMapBadness', aln.badness])
    df = pd.DataFrame(r, columns=['AlignmentId', 'TranscriptId', 'classifier', 'value'])
    df.set_index(['AlignmentId', 'TranscriptId', 'classifier'], inplace=True)
    return df


def paralogy(psl_dict):
    """
    Count the number of occurrences of each parental annotation in the target genome
    :param psl_dict: PslDict from psl module of transMap alignments
    :return: collections.Counter
    """
    r = collections.Counter()
    for aln_id in psl_dict:
        r[tools.nameConversions.strip_alignment_numbers(aln_id)] += 1
    return r


def aln_extends_off_contig(aln):
    """
    Does the alignment extend off of a contig or scaffold?
    aligned: #  unaligned: -  whatever: .  edge: |
             query  |---#####....
             target    |#####....
    OR
    aligned: #  unaligned: -  whatever: .  edge: |
             query  ...######---|
             target ...######|

    :param aln: PslRow object
    :return: boolean
    """
    if aln.t_start == 0 and aln.q_start != 0 or aln.t_end == aln.t_size and aln.q_end != aln.q_size:
        return True
    else:
        return False


def alignment_partial_map(aln):
    """
    Does the query sequence not map entirely?

    a.q_size != a.q_end - a.q_start

    :param aln: PslRow object
    :return: boolean
    """
    return True if aln.q_size != aln.q_end - aln.q_start else False


def aln_abuts_unknown_bases(tx, fasta):
    """
    Do any exons in this alignment immediately touch Ns?

    :param tx: a GenePredTranscript object
    :param fasta: pyfasta Fasta object for genome
    :return: boolean
    """
    chrom = tx.chromosome
    for exon in tx.exon_intervals:
        if exon.start == 0:  # we are at the edge of the contig
            left_base = None
        else:
            left_base = fasta[chrom][exon.start - 1]
        if exon.stop >= len(fasta[chrom]):  # we are at the edge of the contig
            right_base = None
        else:
            right_base = fasta[chrom][exon.stop]
        if left_base == 'N' or right_base == 'N':
            return True
    return False


def aln_contains_unknown_bases(tx, fasta):
    """
    Does this alignment contain unknown bases (Ns)?

    :param tx: a GenePredTranscript object
    :param fasta: pyfasta Fasta object for genome
    :return: boolean
    """
    return 'N' not in tx.get_mrna(fasta)


def long_transcript(tx):
    """
    Is this transcript greater in genomic length than long_tx_size?

    :param tx: a GenePredTranscript object
    :return: boolean
    """
    return True if tx.start - tx.stop >= long_tx_size else False


def synteny(ref_gp_dict, gp_dict):
    """
    Attempts to evaluate the synteny of these transcripts. For each transcript, compares the 5 genes up and down stream
    in the reference genome and counts how many match the transMap results.
    :param ref_gp_dict: Dictionary of GenePredTranscript objects from the reference annotation
    :param gp_dict: Dictionary of GenePredTranscript objects from the transMap output
    :return:
    """
    def create_interval_dict(tx_dict):
        """
        Creates a dict mapping chromosome sequences to gene intervals [chrom][gene_id]: [list of tx intervals]
        Skips huge intervals to avoid mapping issues
        """
        interval_dict = collections.defaultdict(lambda: collections.defaultdict(list))
        for tx in tx_dict.itervalues():
            if len(tx.interval) < long_tx_size:
                interval_dict[tx.chromosome][tx.name2].append(tx.interval)
        return interval_dict

    def merge_interval_dict(interval_dict):
        """Merges the above intervals into the one genic interval."""
        merged_interval_dict = collections.defaultdict(dict)
        for chrom in interval_dict:
            for gene_id, gene_intervals in interval_dict[chrom].iteritems():
                merged_intervals = tools.intervals.gap_merge_intervals(gene_intervals, float('inf'))
                assert len(merged_intervals) == 1
                merged_interval = merged_intervals[0]
                if len(merged_interval) >= long_tx_size:
                    continue
                merged_interval.data = gene_id
                merged_interval_dict[chrom][gene_id] = merged_interval
        return merged_interval_dict

    def sort_interval_dict(merged_interval_dict):
        """Sorts the dict produced by create_interval_dict so that we can do list bisection"""
        sorted_interval_dict = {}
        for chrom in merged_interval_dict:
            sorted_interval_dict[chrom] = sorted(merged_interval_dict[chrom].itervalues())
        return sorted_interval_dict

    def make_ref_interval_map(ref_intervals):
        """Creates a dictionary mapping reference intervals to their name"""
        ref_interval_map = {}
        for interval_list in ref_intervals.itervalues():
            for interval in interval_list:
                assert interval.data not in ref_interval_map
                ref_interval_map[interval.data] = interval
        return ref_interval_map

    # create dictionaries mapping chromosome names to all genic intervals present on the chromosome
    tm_chrom_intervals = sort_interval_dict(merge_interval_dict(create_interval_dict(gp_dict)))
    ref_chrom_intervals = sort_interval_dict(merge_interval_dict(create_interval_dict(ref_gp_dict)))

    # convert the reference to a map that is per-name so that we know where to look
    ref_interval_map = make_ref_interval_map(ref_chrom_intervals)

    # synteny score algorithm
    scores = {}
    for tx in gp_dict.itervalues():
        # find the genes from -5 to +5 in the target genome
        target_intervals = tm_chrom_intervals[tx.chromosome]
        target_position = bisect.bisect_left(target_intervals, tx.interval)
        target_genes = {x.data for x in target_intervals[target_position - 5: target_position + 5]}
        # find the same gene list in the reference genome
        ref_interval = ref_interval_map[tx.name2]
        ref_intervals = ref_chrom_intervals[ref_interval.chromosome]
        ref_position = bisect.bisect_left(ref_intervals, ref_interval)
        reference_genes = {x.data for x in ref_intervals[ref_position - 5: ref_position + 5]}
        scores[tx.name] = len(reference_genes & target_genes)
    return scores
