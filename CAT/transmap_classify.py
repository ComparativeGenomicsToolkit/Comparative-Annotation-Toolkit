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
8. Distance: What is the tree building distance from the parent gene?
"""
import bisect
import collections
import argparse
import ete3
import itertools

import pandas as pd
from toil.common import Toil
from toil.job import Job

import tools.bio
import tools.nameConversions
import tools.psl
import tools.dataOps
import tools.fileOps
import tools.transcripts
import tools.toilInterface
import tools.procOps

# hard coded variables
# hard coded long transMap size. Bigger than 3 megabases is probably a spurious alignment.
long_tx_size = 3 * 10 ** 6


def transmap_classify(tm_eval_args, toil_options):
    """
    Runs alignment classification based on transMap PSLs, genePreds and the genome FASTA. Launches a toil pipeline
    within this module which runs MUSCLE + FastTree on all alignments, reporting the genetic distance
    :param tm_eval_args: argparse Namespace produced by EvaluateTransMap.get_args()
    :param toil_options: toil options Namespace object
    :return: DataFrame
    """
    psl_dict = tools.psl.get_alignment_dict(tm_eval_args.tm_psl)
    gp_dict = tools.transcripts.get_gene_pred_dict(tm_eval_args.tm_gp)
    ref_gp_dict = tools.transcripts.get_gene_pred_dict(tm_eval_args.annotation_gp)
    fasta = tools.bio.get_sequence_dict(tm_eval_args.fasta)

    paralog_count, paralog_names = paralogy(psl_dict)  # we have to count paralogs globally
    distances = distance_toil_pipeline(tm_eval_args, toil_options, paralog_names)

    synteny_scores = synteny(ref_gp_dict, gp_dict)  # we also have to score synteny globally

    r = []
    for aln_id, tx in gp_dict.iteritems():
        aln = psl_dict[aln_id]
        tx_id = tools.nameConversions.strip_alignment_numbers(aln_id)
        gene_id = ref_gp_dict[tx_id].name2
        r.append([aln_id, tx_id, gene_id, 'Paralogy', paralog_count[tools.nameConversions.strip_alignment_numbers(aln_id)]])
        r.append([aln_id, tx_id, gene_id, 'Distance', distances[aln_id]])
        r.append([aln_id, tx_id, gene_id, 'Synteny', synteny_scores[aln_id]])
        r.append([aln_id, tx_id, gene_id, 'AlnExtendsOffContig', aln_extends_off_contig(aln)])
        r.append([aln_id, tx_id, gene_id, 'AlnPartialMap', alignment_partial_map(aln)])
        r.append([aln_id, tx_id, gene_id, 'AlnAbutsUnknownBases', aln_abuts_unknown_bases(tx, fasta)])
        r.append([aln_id, tx_id, gene_id, 'AlnContainsUnknownBases', aln_contains_unknown_bases(tx, fasta)])
        r.append([aln_id, tx_id, gene_id, 'LongTranscript', long_transcript(tx)])
        r.append([aln_id, tx_id, gene_id, 'TransMapCoverage', aln.coverage])
        r.append([aln_id, tx_id, gene_id, 'TransMapIdentity', aln.identity])
        r.append([aln_id, tx_id, gene_id, 'TransMapBadness', aln.badness])
    df = pd.DataFrame(r, columns=['AlignmentId', 'TranscriptId', 'GeneId', 'classifier', 'value'])
    df.value = pd.to_numeric(df.value)
    df.set_index(['AlignmentId', 'TranscriptId', 'GeneId', 'classifier'], inplace=True)
    return df


###
# Classifiers
###


def paralogy(psl_dict):
    """
    Count the number of occurrences of each parental annotation in the target genome
    :param psl_dict: PslDict from psl module of transMap alignments
    :return: collections.Counter, collections.defaultdict
    """
    counts = collections.Counter()
    names = collections.defaultdict(list)
    for aln_id in psl_dict:
        counts[tools.nameConversions.strip_alignment_numbers(aln_id)] += 1
        names[tools.nameConversions.strip_alignment_numbers(aln_id)].append(aln_id)
    return counts, names


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
    return 'N' in tx.get_mrna(fasta)


def long_transcript(tx):
    """
    Is this transcript greater in genomic length than long_tx_size?

    :param tx: a GenePredTranscript object
    :return: boolean
    """
    return True if tx.start - tx.stop >= long_tx_size else False


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
        # find the genes from -3 to +3 in the target genome
        target_intervals = tm_chrom_intervals[tx.chromosome]
        target_position = bisect.bisect_left(target_intervals, tx.interval)
        target_genes = {x.data for x in target_intervals[target_position - 3: target_position + 3]}
        # find the same gene list in the reference genome
        ref_interval = ref_interval_map[tx.name2]
        ref_intervals = ref_chrom_intervals[ref_interval.chromosome]
        ref_position = bisect.bisect_left(ref_intervals, ref_interval)
        reference_genes = {x.data for x in ref_intervals[ref_position - 3: ref_position + 3]}
        scores[tx.name] = len(reference_genes & target_genes)
    return scores


###
# Paralog count and distance measuring toil pipeline
###

def distance_toil_pipeline(args, toil_options, paralog_dict):
    """
    Produces a dictionary of
    :param args: argparse Namespace produced by EvaluateTransMap.get_args()
    :param toil_options: toil options Namespace object
    :param paralog_dict: dictionary mapping paralogous transcripts
    :return: dictionary mapping aln_ids to measured distances
    """
    with Toil(toil_options) as toil:
        if not toil.options.restart:
            input_file_ids = argparse.Namespace()
            input_file_ids.ref_genome_fasta = tools.toilInterface.write_fasta_to_filestore(toil, args.ref_fasta)
            input_file_ids.genome_fasta = tools.toilInterface.write_fasta_to_filestore(toil, args.fasta)
            input_file_ids.annotation_gp = toil.importFile('file://' + args.annotation_gp)
            input_file_ids.tm_gp = toil.importFile('file://' + args.tm_gp)
            job = Job.wrapJobFn(setup, paralog_dict, input_file_ids, memory='16G')
            distances = toil.start(job)
        else:
            distances = toil.restart()
        return distances


def setup(job, paralog_dict, input_file_ids):
    # import the files
    tm_gp = job.fileStore.readGlobalFile(input_file_ids.tm_gp)
    annotation_gp = job.fileStore.readGlobalFile(input_file_ids.annotation_gp)
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.genome_fasta,
                                                                 prefix='genome', upper=False)
    ref_genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.ref_genome_fasta,
                                                                     prefix='ref_genome', upper=False)
    tx_dict = tools.transcripts.get_gene_pred_dict(tm_gp)
    ref_tx_dict = tools.transcripts.get_gene_pred_dict(annotation_gp)
    # dists will hold a list of lists of results to be collapsed in next job
    distances = []
    tx_iter = aln_seq_iter(paralog_dict, tx_dict, ref_tx_dict, genome_fasta, ref_genome_fasta)
    for chunk, mem in group_transcripts(tx_iter):
        distances.append(job.addChildJobFn(calculate_distances, chunk, memory=mem).rv())
    return job.addFollowOnJobFn(merge, distances).rv()


def calculate_distances(job, chunk):
    """
    Calculates the phylogenetic distance for this chunk
    :param chunk: chunk of (tx_id, tx_list, fasta_str) tuples
    :return: list of [aln_id: distance] lists
    """
    tmp_fasta = tools.fileOps.get_tmp_toil_file(suffix='fa')
    tmp_aligned_fasta = tools.fileOps.get_tmp_toil_file(suffix='aligned.fa')
    distances = []
    for tx_id, seq_list in chunk:
        with open(tmp_fasta, 'w') as outf:
            for aln_id, seq in seq_list:
                tools.bio.write_fasta(outf, aln_id, seq)
        cmd = ['muscle', '-quiet', '-in', tmp_fasta, '-diags1', '-sv', '-maxiters', '2', '-out', tmp_aligned_fasta]
        tools.procOps.run_proc(cmd)
        cmd = ['FastTree', '-nt', '-gtr', '-quiet', '-pseudo', tmp_aligned_fasta]
        newick = tools.procOps.call_proc_lines(cmd)[0]
        distances.extend(parse_tree(newick, tx_id))
    return distances


def merge(job, distances):
    """
    Merges the final distances, reformats paralogs
    :param distances: list of lists of distance measures
    :return: distance_dict
    """
    # flatten distances into dictionary
    distance_dict = {}
    for l in distances:
        for aln_id, d in l:
            distance_dict[aln_id] = d
    return distance_dict


###
# Helper functions for toil pipeline
###


def aln_seq_iter(paralog_dict, tx_dict, ref_tx_dict, genome_fasta, ref_genome_fasta):
    """yields tuples of (tx_id, [aln_ids], fasta_string for alignment and tree building"""
    for tx_id, tx_list in paralog_dict.iteritems():
        ref_tx = ref_tx_dict[tx_id]
        seq_list = [[tx_id, ref_tx.get_mrna(ref_genome_fasta)]]
        for aln_id in tx_list:
            tx = tx_dict[aln_id]
            seq_list.append([aln_id, tx.get_mrna(genome_fasta)])
        yield tx_id, seq_list


def group_transcripts(tx_iter, num_bases=0.20 * 10 ** 6, max_seqs=100):
    """
    Group up transcripts by num_bases, unless that exceeds max_seqs. A greedy implementation of the bin packing problem.

    Memory requirements are calculated by the largest pair of sequences. num_bases is by the largest sequence.
    """
    tx_id, seq_list = tx_iter.next()
    this_bin = [(tx_id, seq_list)]
    num_seqs = 1
    bin_base_count = max(len(seq) for name, seq in seq_list)
    mem = find_memory_requirements(seq_list)
    for tx_id, seq_list in tx_iter:
        bin_base_count += max(len(seq) for name, seq in seq_list)
        num_seqs += 1
        mem = max(mem, find_memory_requirements(seq_list))
        if bin_base_count >= num_bases or num_seqs >= max_seqs:
            yield this_bin, '{}G'.format(mem)
            this_bin = [(tx_id, seq_list)]
            bin_base_count = max(len(seq) for name, seq in seq_list)
            num_seqs = 1
            mem = find_memory_requirements(seq_list)
        else:
            this_bin.append((tx_id, seq_list))
    yield this_bin, '{}G'.format(mem)


def find_memory_requirements(seq_list):
    """
    Determines how much memory a pair of transcripts will need for alignment. Does this by multplying their sizes,
    rounding up to the nearest multiple of 4 with 2 extra gigabytes of headroom
    """
    def calc_diag(s1, s2):
        return 2 + 1.0 * len(s1) * len(s2) / 10 ** 9
    seqs = [seq for name, seq in seq_list]
    requirement = max(calc_diag(*x) for x in itertools.combinations(seqs, 2))
    mem = requirement if requirement % 4 == 0 else requirement + 4 - requirement % 4
    return int(mem)


def parse_tree(newick, tx_id):
    """parses the newick tree from FastTree"""
    t = ete3.Tree(newick)
    return [[node.name, t.get_distance(tx_id, node)] for node in t]


