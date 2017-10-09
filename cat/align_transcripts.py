"""
Toil pipeline to align all transcripts in a source genePred to a target genePred.

Alignment is performed two ways -- full mRNA and in-frame CDS. Alignment is done for all TM/TMR transcripts and all
protein coding transMap transcripts.
"""
import argparse
import collections
import itertools
import logging

from toil.fileStore import FileID
from toil.common import Toil
from toil.job import Job

import tools.bio
import tools.dataOps
import tools.fileOps
import tools.nameConversions
import tools.pipeline
import tools.procOps
import tools.psl
import tools.sqlInterface
import tools.toilInterface
import tools.transcripts


def align_transcripts(args, toil_options):
    """
    Main entry function for transcript alignment toil pipeline
    :param args: dictionary of arguments from CAT
    :param toil_options: toil options Namespace object
    """
    with Toil(toil_options) as t:
        if not t.options.restart:
            input_file_ids = argparse.Namespace()
            input_file_ids.ref_genome_fasta = tools.toilInterface.write_fasta_to_filestore(t, args.ref_genome_fasta)
            input_file_ids.genome_fasta = tools.toilInterface.write_fasta_to_filestore(t, args.genome_fasta)
            input_file_ids.annotation_gp = FileID.forPath(t.importFile('file://' + args.annotation_gp),
                                                          args.annotation_gp)
            input_file_ids.ref_db = FileID.forPath(t.importFile('file://' + args.ref_db_path), args.ref_db_path)
            input_file_ids.modes = {}
            file_ids = [input_file_ids.ref_genome_fasta, input_file_ids.genome_fasta, input_file_ids.annotation_gp,
                        input_file_ids.ref_db]
            for mode in args.transcript_modes:
                input_file_ids.modes[mode] = t.importFile('file://' + args.transcript_modes[mode]['gp'])
                file_ids.append(input_file_ids.modes[mode])
            disk_usage = tools.toilInterface.find_total_disk_usage(file_ids)
            job = Job.wrapJobFn(setup, args, input_file_ids, memory='16G', disk=disk_usage)
            results_file_ids = t.start(job)
        else:
            results_file_ids = t.restart()
        for file_path, file_id in results_file_ids.iteritems():
            tools.fileOps.ensure_file_dir(file_path)
            t.exportFile(file_id, 'file://' + file_path)


def setup(job, args, input_file_ids):
    """
    First function for align_transcripts pipeline. Splits up the genePred entries into chunks that will be aligned
    with BLAT.
    :param args: dictionary of arguments from CAT
    :param input_file_ids: dictionary of fileStore file IDs for the inputs to this pipeline
    """
    job.fileStore.logToMaster('Beginning Align Transcripts run on {}'.format(args.genome), level=logging.INFO)
    # load all fileStore files necessary
    annotation_gp = job.fileStore.readGlobalFile(input_file_ids.annotation_gp)
    ref_genome_db = job.fileStore.readGlobalFile(input_file_ids.ref_db)
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.genome_fasta,
                                                                 prefix='genome', upper=False)
    ref_genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.ref_genome_fasta,
                                                                     prefix='ref_genome', upper=False)
    # load required reference data into memory
    tx_biotype_map = tools.sqlInterface.get_transcript_biotype_map(ref_genome_db)
    ref_transcript_dict = tools.transcripts.get_gene_pred_dict(annotation_gp)
    # will hold a mapping of output file paths to lists of Promise objects containing output
    results = collections.defaultdict(list)
    for tx_mode in ['transMap', 'augTM', 'augTMR']:
        if tx_mode not in args.transcript_modes:
            continue
        # output file paths
        mrna_path = args.transcript_modes[tx_mode]['mRNA']
        cds_path = args.transcript_modes[tx_mode]['CDS']
        # begin loading transcripts and sequences
        gp_path = job.fileStore.readGlobalFile(input_file_ids.modes[tx_mode])
        transcript_dict = tools.transcripts.get_gene_pred_dict(gp_path)
        transcript_dict = {aln_id: tx for aln_id, tx in transcript_dict.iteritems() if
                           tx_biotype_map[tools.nameConversions.strip_alignment_numbers(aln_id)] == 'protein_coding'}
        for aln_mode, out_path in zip(*[['mRNA', 'CDS'], [mrna_path, cds_path]]):
            seq_iter = get_alignment_sequences(transcript_dict, ref_transcript_dict, genome_fasta,
                                               ref_genome_fasta, aln_mode)
            for chunk in group_transcripts(seq_iter):
                j = job.addChildJobFn(run_blat_chunk, chunk, aln_mode, memory='8G', disk='2G')
                results[out_path].append(j.rv())

    if len(results) == 0:
        err_msg = 'Align Transcripts pipeline did not detect any input genePreds for {}'.format(args.genome)
        raise RuntimeError(err_msg)
    # convert the results Promises into resolved values
    return job.addFollowOnJobFn(merge, results, args, memory='8G', disk='4G').rv()


def get_alignment_sequences(transcript_dict, ref_transcript_dict, genome_fasta, ref_genome_fasta, mode):
    """Generator that yields a tuple of (tx_id, tx_seq, ref_tx_id, ref_tx_seq)"""
    assert mode in ['mRNA', 'CDS']
    for tx_id, tx in transcript_dict.iteritems():
        ref_tx_id = tools.nameConversions.strip_alignment_numbers(tx_id)
        ref_tx = ref_transcript_dict[ref_tx_id]
        tx_seq = tx.get_mrna(genome_fasta) if mode == 'mRNA' else tx.get_cds(genome_fasta)
        ref_tx_seq = ref_tx.get_mrna(ref_genome_fasta) if mode == 'mRNA' else ref_tx.get_cds(ref_genome_fasta)
        if len(ref_tx_seq) > 50 and len(tx_seq) > 50:
            yield tx_id, tx_seq, ref_tx_id, ref_tx_seq


def run_blat_chunk(job, chunk, mode):
    """
    Runs an alignment chunk through BLAT for either coding or non-coding transcripts
    :param chunk: List of (tx_id, tx_seq, ref_tx_id, ref_tx_seq) tuples
    :param mode: One of ['mRNA', 'CDS']. Determines what mode of alignment we will perform.
    :return: List of PSL output
    """
    def parse_blat(tmp_psl):
        # filter for only + alignments, as we are expecting to be on the same strand
        # translation alignments have explicit strand, and we only want ++
        filter_strand = '+' if mode == 'mRNA' else '++'
        psls = [psl for psl in tools.psl.psl_iterator(tmp_psl) if psl.strand == filter_strand]
        if len(psls) == 0:
            return None
        longest = sorted(psls, key=lambda p: -p.coverage)[0]
        return '\t'.join(longest.psl_string())

    assert mode in ['mRNA', 'CDS']
    tmp_ref = tools.fileOps.get_tmp_toil_file()
    tmp_tgt = tools.fileOps.get_tmp_toil_file()
    tmp_psl = tools.fileOps.get_tmp_toil_file()
    tmp_filtered_psl = tools.fileOps.get_tmp_toil_file()
    results = []
    if mode == 'mRNA':
        cmd = ['blat', '-noHead', '-minIdentity=0', tmp_ref, tmp_tgt, tmp_psl]
    else:  # mode == CDS. Filter these for problematic alignments that happen in edge cases
        cmd = ['blat', '-t=dnax', '-q=rnax', '-noHead', '-minIdentity=0', tmp_ref, tmp_tgt, tmp_psl]
    for tx_id, tx_seq, ref_tx_id, ref_tx_seq in chunk:
        with open(tmp_ref, 'w') as tmp_ref_h:
            tools.bio.write_fasta(tmp_ref_h, ref_tx_id, ref_tx_seq)
        with open(tmp_tgt, 'w') as tmp_tgt_h:
            tools.bio.write_fasta(tmp_tgt_h, tx_id, tx_seq)
        tools.procOps.run_proc(cmd)
        try:
            tools.procOps.run_proc(['pslCheck', '-quiet', tmp_psl, '-pass={}'.format(tmp_filtered_psl)])
        except tools.pipeline.ProcException:
            pass
        results.append(parse_blat(tmp_filtered_psl))
    return results


def merge(job, results, args):
    """
    Merge together chain files.
    :param results: dict of list of promises from each alignment chunk for each category
    :param args: arguments to the pipeline
    :return:
    """
    job.fileStore.logToMaster('Merging Alignment output for {}'.format(args.genome), level=logging.INFO)
    results_file_ids = {}
    for gp_category, result_list in results.iteritems():
        tmp_results_file = tools.fileOps.get_tmp_toil_file()
        with open(tmp_results_file, 'w') as outf:
            for line in itertools.chain.from_iterable(result_list):  # results is list of lists
                if line is not None:
                    outf.write(line + '\n')
        results_file_ids[gp_category] = job.fileStore.writeGlobalFile(tmp_results_file)
    return results_file_ids


###
# Helper functions
###


def group_transcripts(tx_iter, num_bases=10 ** 6, max_seqs=1000):
    """
    Group up transcripts by num_bases, unless that exceeds max_seqs. A greedy implementation of the bin packing problem.
    Helps speed up the execution of BLAT when faced with very large genes
    """
    tx_id, tx_seq, ref_tx_id, ref_tx_seq = tx_iter.next()
    this_bin = [(tx_id, tx_seq, ref_tx_id, ref_tx_seq)]
    bin_base_count = len(tx_seq)
    num_seqs = 1
    for tx_id, tx_seq, ref_tx_id, ref_tx_seq in tx_iter:
        bin_base_count += len(tx_seq)
        num_seqs += 1
        if bin_base_count >= num_bases or num_seqs >= max_seqs:
            yield this_bin
            this_bin = [(tx_id, tx_seq, ref_tx_id, ref_tx_seq)]
            bin_base_count = len(tx_seq)
            num_seqs = 1
        else:
            this_bin.append((tx_id, tx_seq, ref_tx_id, ref_tx_seq))
    yield this_bin
