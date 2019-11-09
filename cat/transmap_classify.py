"""
Classify transMap transcripts producing the TransMapEvaluation table for each genome's database

1. AlnExtendsOffConfig: Does this alignment run off the end of a contig?
2. AlnPartialMap: Did this transcript not map completely?
3. AlnAbutsUnknownBases: Does this alignment have Ns immediately touching any exons?
4. PercentN: Percent of bases aligned to Ns
5. TransMapCoverage
6. TransMapIdentity
7. TransMapGoodness
8. TransMapOriginalIntronsPercent: The number of transMap introns within a wiggle distance of a intron in the parent
 transcript in transcript coordinates.
9. Synteny. Count of the # of genes that match the reference in both directions (+/- 5 genes)
10. ValidStart -- start with ATG?
11. ValidStop -- valid stop codon (in frame)?
12. ProperOrf -- is the orf a multiple of 3?
"""
import bisect
import collections
import pandas as pd

import tools.bio
import tools.nameConversions
import tools.psl
import tools.dataOps
import tools.fileOps
import tools.transcripts
import tools.toilInterface
import tools.procOps
import tools.tm2hints
import tools.mathOps


def transmap_classify(tm_eval_args):
    """
    Wrapper function that runs alignment classification based on transMap PSLs, genePreds and the genome FASTA.
    :param tm_eval_args: argparse Namespace produced by EvaluateTransMap.get_args()
    :return: DataFrame
    """
    psl_dict = tools.psl.get_alignment_dict(tm_eval_args.filtered_tm_psl)
    ref_psl_dict = tools.psl.get_alignment_dict(tm_eval_args.ref_psl)
    gp_dict = tools.transcripts.get_gene_pred_dict(tm_eval_args.filtered_tm_gp)
    ref_gp_dict = tools.transcripts.get_gene_pred_dict(tm_eval_args.annotation_gp)
    fasta = tools.bio.get_sequence_dict(tm_eval_args.fasta)

    synteny_scores = synteny(ref_gp_dict, gp_dict)

    r = []
    for aln_id, tx in gp_dict.items():
        aln = psl_dict[aln_id]
        tx_id = tools.nameConversions.strip_alignment_numbers(aln_id)
        ref_aln = ref_psl_dict[tx_id]
        gene_id = ref_gp_dict[tx_id].name2
        r.append([aln_id, tx_id, gene_id, 'AlnExtendsOffContig', aln_extends_off_contig(aln)])
        r.append([aln_id, tx_id, gene_id, 'AlnPartialMap', alignment_partial_map(aln)])
        r.append([aln_id, tx_id, gene_id, 'AlnAbutsUnknownBases', aln_abuts_unknown_bases(tx, fasta)])
        r.append([aln_id, tx_id, gene_id, 'PercentN', aln.percent_n])
        r.append([aln_id, tx_id, gene_id, 'TransMapCoverage', 100 * aln.coverage])
        r.append([aln_id, tx_id, gene_id, 'TransMapIdentity', 100 * aln.identity])
        r.append([aln_id, tx_id, gene_id, 'TransMapGoodness', 100 * (1 - aln.badness)])
        r.append([aln_id, tx_id, gene_id, 'TransMapOriginalIntronsPercent', percent_original_introns(aln, tx, ref_aln)])
        r.append([aln_id, tx_id, gene_id, 'Synteny', synteny_scores[aln_id]])
        r.append([aln_id, tx_id, gene_id, 'ValidStart', tools.transcripts.has_start_codon(fasta, tx)])
        r.append([aln_id, tx_id, gene_id, 'ValidStop', tools.transcripts.has_stop_codon(fasta, tx)])
        r.append([aln_id, tx_id, gene_id, 'ProperOrf', tx.cds_size % 3 == 0])
    df = pd.DataFrame(r, columns=['AlignmentId', 'TranscriptId', 'GeneId', 'classifier', 'value'])
    df.value = pd.to_numeric(df.value)
    return df.set_index(['GeneId', 'TranscriptId', 'AlignmentId', 'classifier'])


###
# Classifiers
###


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
    return aln.q_size != aln.q_end - aln.q_start


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


def synteny(ref_gp_dict, gp_dict):
    """
    Attempts to evaluate the synteny of these transcripts. For each transcript, compares the 3 genes up and down stream
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
        for tx in tx_dict.values():
            interval_dict[tx.chromosome][tx.name2].append(tx.interval)
        return interval_dict

    def merge_interval_dict(interval_dict):
        """Merges the above intervals into the one genic interval."""
        merged_interval_dict = collections.defaultdict(dict)
        for chrom in interval_dict:
            for gene_id, gene_intervals in interval_dict[chrom].items():
                merged_intervals = tools.intervals.gap_merge_intervals(gene_intervals, float('inf'))
                assert len(merged_intervals) == 1
                merged_interval = merged_intervals[0]
                merged_interval.data = gene_id
                merged_interval_dict[chrom][gene_id] = merged_interval
        return merged_interval_dict

    def sort_interval_dict(merged_interval_dict):
        """Sorts the dict produced by create_interval_dict so that we can do list bisection"""
        sorted_interval_dict = {}
        for chrom in merged_interval_dict:
            sorted_interval_dict[chrom] = sorted(merged_interval_dict[chrom].values())
        return sorted_interval_dict

    def make_ref_interval_map(ref_intervals):
        """Creates a dictionary mapping reference intervals to their name"""
        ref_interval_map = {}
        for interval_list in ref_intervals.values():
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
    for tx in gp_dict.values():
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


def percent_original_introns(aln, tx, ref_aln):
    """
    Calculates the intron support vector, using code from tm2hints, but shrinking the fuzz distance to match the
    alignment classifiers.
    Returns the number of introns that are within wiggle distance
    :param aln: PslRow object representing the transMapped transcript
    :param tx: GenePredTranscript object representing the transMapped transcript
    :param ref_aln: PslRow object representing the reference transcript
    :return: float between 0 and 100
    """
    ref_starts = tools.tm2hints.fix_ref_q_starts(ref_aln)
    c = 0
    for i in tx.intron_intervals:
        if tools.tm2hints.is_fuzzy_intron(i, aln, ref_starts, fuzz_distance=7):
            c += 1
    return 100 * tools.mathOps.format_ratio(c, len(tx.intron_intervals), resolve_nan=None)
