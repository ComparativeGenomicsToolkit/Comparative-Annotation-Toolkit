"""
Runs AugustusPB on a target genome
"""

import argparse
import collections
from toil.common import Toil
from toil.job import Job

import tools.bio
import tools.misc
import tools.dataOps
import tools.fileOps
import tools.intervals
import tools.procOps
import tools.sqlInterface
import tools.transcripts
import tools.parentGeneAssignment
import tools.toilInterface


###
# AugustusPB pipeline section
###


def augustus_pb(args, toil_options):
    """
    Main entry function for AugustusPB toil pipeline
    :param args: dictionary of arguments from CAT
    :param toil_options: toil options Namespace object
    :return:
    """
    with Toil(toil_options) as toil:
        if not toil.options.restart:
            input_file_ids = argparse.Namespace()
            input_file_ids.genome_fasta = tools.toilInterface.write_fasta_to_filestore(toil, args.genome_fasta)
            input_file_ids.filtered_tm_gp = toil.importFile('file://' + args.filtered_tm_gp)
            input_file_ids.unfiltered_tm_gp = toil.importFile('file://' + args.unfiltered_tm_gp)
            input_file_ids.ref_db_path = toil.importFile('file://' + args.ref_db_path)
            input_file_ids.pb_cfg = toil.importFile('file://' + args.pb_cfg)
            input_file_ids.hints_gff = toil.importFile('file://' + args.hints_gff)
            job = Job.wrapJobFn(setup, args, input_file_ids, memory='8G')
            raw_gtf_file_id, (gtf_file_id, df, fail_count) = toil.start(job)
        else:
            raw_gtf_file_id, (gtf_file_id, df, fail_count) = toil.restart()
        tools.fileOps.ensure_file_dir(args.augustus_pb_raw_gtf)
        toil.exportFile(raw_gtf_file_id, 'file://' + args.augustus_pb_raw_gtf)
        tools.fileOps.ensure_file_dir(args.augustus_pb_gtf)
        toil.exportFile(gtf_file_id, 'file://' + args.augustus_pb_gtf)
        return df, fail_count


def setup(job, args, input_file_ids):
    """
    Entry function for running AugustusPB.
    The genome is chunked up and the resulting gene sets merged using joingenes.
    """
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.genome_fasta,
                                                                 prefix='genome', upper=False)

    # calculate overlapping intervals. If the final interval is small (<= 50% of total interval size), merge it
    intervals = collections.defaultdict(list)
    for chrom in genome_fasta:
        chrom_size = len(genome_fasta[chrom])
        for start in xrange(0, chrom_size, args.chunksize - args.overlap):
            stop = min(start + args.chunksize, chrom_size)
            intervals[chrom].append([start, stop])

    for chrom, interval_list in intervals.iteritems():
        if len(interval_list) < 2:
            continue
        last_start, last_stop = interval_list[-1]
        if last_stop - last_start <= 0.5 * args.chunksize:
            del interval_list[-1]
            interval_list[-1][-1] = last_stop

    predictions = []
    for chrom, interval_list in intervals.iteritems():
        for start, stop in interval_list:
            j = job.addChildJobFn(augustus_pb_chunk, args, input_file_ids, chrom, start, stop)
            predictions.append(j.rv())

    # results contains a 3 member tuple of [gff_file_id, dataframe, fail_count]
    # where the dataframe contains the alternative parental txs and fail_count is the # of transcripts discarded
    results = job.addFollowOnJobFn(join_genes, input_file_ids, predictions, memory='8G').rv()
    return results


def augustus_pb_chunk(job, args, input_file_ids, chrom, start, stop):
    """
    core function that runs AugustusPB on one genome chunk
    """
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.genome_fasta,
                                                                 prefix='genome', upper=False)
    hints = job.fileStore.readGlobalFile(input_file_ids.hints_gff)

    # slice out only the relevant hints
    hints_subset = tools.fileOps.get_tmp_toil_file()
    cmd = ['awk', '($1 == "{}" && $4 >= {} && $5 <= {}) {{print $0}}'.format(chrom, start, stop)]
    tools.procOps.run_proc(cmd, stdin=hints, stdout=hints_subset)

    pb_cfg = job.fileStore.readGlobalFile(input_file_ids.pb_cfg)
    tmp_fasta = tools.fileOps.get_tmp_toil_file()
    tools.bio.write_fasta(tmp_fasta, chrom, genome_fasta[chrom][start:stop])
    results = tools.fileOps.get_tmp_toil_file()

    cmd = ['augustus', '--UTR=1', '--softmasking=1', '--allow_hinted_splicesites=atac',
           '--alternatives-from-evidence=1',
           '--hintsfile={}'.format(hints_subset),
           '--extrinsicCfgFile={}'.format(pb_cfg),
           '--species={}'.format(args.species),
           '--predictionStart=-{}'.format(start), '--predictionEnd=-{}'.format(stop),
           tmp_fasta]
    tools.procOps.run_proc(cmd, stdout=results)
    return job.fileStore.writeGlobalFile(results)


def join_genes(job, input_file_ids, predictions):
    """
    uses the auxiliary tool 'joingenes' from the
    Augustus package to intelligently merge gene sets
    - removes duplicated Txs or truncated Txs that are contained in other Txs (trivial)
    - fixes truncated Txs at alignment boundaries,
      e.g. by merging them with other Txs (non trivial, introduces new Txs)

    Calls out to the parental gene assignment pipeline
    """
    raw_gtf_file = tools.fileOps.get_tmp_toil_file()
    raw_gtf_fofn = tools.fileOps.get_tmp_toil_file()
    with open(raw_gtf_file, 'w') as raw_handle, open(raw_gtf_fofn, 'w') as fofn_handle:
        for chunk in predictions:
            local_path = job.fileStore.readGlobalFile(chunk)
            fofn_handle.write(local_path + '\n')
            for line in open(local_path):
                raw_handle.write(line)

    join_genes_file = tools.fileOps.get_tmp_toil_file()
    cmd = [['joingenes', '-f', raw_gtf_fofn, '-o', '/dev/stdout', '--alternatives'],
           ['grep', '-P', '\tAUGUSTUS\t(exon|CDS|start_codon|stop_codon|tts|tss)\t'],
           ['sed', ' s/jg/augPB-/g']]
    tools.procOps.run_proc(cmd, stdout=join_genes_file)
    joined_file_id = job.fileStore.writeGlobalFile(join_genes_file)
    raw_gtf_file_id = job.fileStore.writeGlobalFile(raw_gtf_file)
    j = job.addFollowOnJobFn(tools.parentGeneAssignment.assign_parents, input_file_ids.ref_db_path,
                             input_file_ids.filtered_tm_gp, input_file_ids.unfiltered_tm_gp,
                             joined_file_id, 'AugustusPB', memory='8G')
    return raw_gtf_file_id, j.rv()
