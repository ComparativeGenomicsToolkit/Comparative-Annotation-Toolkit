"""
Runs AugustusPB on a target genome
"""
import argparse
import collections
import os

from toil.fileStore import FileID
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
    with Toil(toil_options) as t:
        if not t.options.restart:
            input_file_ids = argparse.Namespace()
            input_file_ids.genome_fasta = tools.toilInterface.write_fasta_to_filestore(t, args.genome_fasta)
            input_file_ids.chrom_sizes = FileID.forPath(t.importFile('file://' + args.chrom_sizes), args.chrom_sizes)
            input_file_ids.pb_cfg = FileID.forPath(t.importFile('file://' + args.pb_cfg), args.pb_cfg)
            input_file_ids.hints_gff = FileID.forPath(t.importFile('file://' + args.hints_gff), args.hints_gff)
            job = Job.wrapJobFn(setup, args, input_file_ids, memory='16G', disk='32G')
            raw_gtf_file_id, gtf_file_id, joined_gp_file_id = t.start(job)
        else:
            raw_gtf_file_id, gtf_file_id, joined_gp_file_id = t.restart()
        tools.fileOps.ensure_file_dir(args.augustus_pb_raw_gtf)
        t.exportFile(raw_gtf_file_id, 'file://' + args.augustus_pb_raw_gtf)
        t.exportFile(gtf_file_id, 'file://' + args.augustus_pb_gtf)
        t.exportFile(joined_gp_file_id, 'file://' + args.augustus_pb_gp)


def setup(job, args, input_file_ids):
    """
    Entry function for running AugustusPB.
    The genome is chunked up and the resulting gene sets merged using joingenes.
    """
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.genome_fasta,
                                                                 prefix='genome', upper=False)

    # load only PB hints
    hints_file = job.fileStore.readGlobalFile(input_file_ids.hints_gff)
    hints = [x.split('\t') for x in open(hints_file) if 'src=PB' in x]

    if len(hints) == 0:
        raise RuntimeError('No PB hints found.')

    # convert the start/stops to ints
    # break up by chromosome
    hints_by_chrom = collections.defaultdict(list)
    for h in hints:
        h[3] = int(h[3])
        h[4] = int(h[4])
        hints_by_chrom[h[0]].append(h)

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
            hints = [h for h in hints_by_chrom[chrom] if h[3] >= start and h[4] <= stop]
            if len(hints) == 0:
                continue  # no reason to compute an empty chunk
            tmp_hints = tools.fileOps.get_tmp_toil_file()
            with open(tmp_hints, 'w') as outf:
                for h in hints:
                    tools.fileOps.print_row(outf, h)
            hints_file_id = job.fileStore.writeGlobalFile(tmp_hints)
            j = job.addChildJobFn(augustus_pb_chunk, args, input_file_ids, hints_file_id, chrom, start, stop,
                                  memory='8G', disk='8G')
            predictions.append(j.rv())

    # results contains a 3 member tuple of [raw_gtf_file_id, gtf_file_id, joined_gp_file_id]
    results = job.addFollowOnJobFn(join_genes, predictions, memory='8G', disk='8G').rv()
    return results


def augustus_pb_chunk(job, args, input_file_ids, hints_file_id, chrom, start, stop):
    """
    core function that runs AugustusPB on one genome chunk
    """
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.genome_fasta,
                                                                 prefix='genome', upper=False)
    hints = job.fileStore.readGlobalFile(hints_file_id)
    pb_cfg = job.fileStore.readGlobalFile(input_file_ids.pb_cfg)
    tmp_fasta = tools.fileOps.get_tmp_toil_file()
    tools.bio.write_fasta(tmp_fasta, chrom, genome_fasta[chrom][start:stop])
    results = tools.fileOps.get_tmp_toil_file()

    cmd = ['augustus', '--softmasking=1', '--allow_hinted_splicesites=atac',
           '--alternatives-from-evidence=1', '--UTR={}'.format(int(args.utr)),
           '--hintsfile={}'.format(hints),
           '--extrinsicCfgFile={}'.format(pb_cfg),
           '--species={}'.format(args.species),
           '--/augustus/verbosity=0',
           '--predictionStart=-{}'.format(start), '--predictionEnd=-{}'.format(start),
           tmp_fasta]
    tools.procOps.run_proc(cmd, stdout=results)
    return job.fileStore.writeGlobalFile(results)


def join_genes(job, gff_chunks):
    """
    uses the auxiliary tool 'joingenes' from the
    Augustus package to intelligently merge gene sets
    - removes duplicated Txs or truncated Txs that are contained in other Txs (trivial)
    - fixes truncated Txs at alignment boundaries,
      e.g. by merging them with other Txs (non trivial, introduces new Txs)
    """
    raw_gtf_file = tools.fileOps.get_tmp_toil_file()
    raw_gtf_fofn = tools.fileOps.get_tmp_toil_file()
    with open(raw_gtf_file, 'w') as raw_handle, open(raw_gtf_fofn, 'w') as fofn_handle:
        for chunk in gff_chunks:
            local_path = job.fileStore.readGlobalFile(chunk)
            if os.environ.get('CAT_BINARY_MODE') == 'singularity':
                local_path = tools.procOps.singularify_arg(local_path)
            fofn_handle.write(local_path + '\n')
            for line in open(local_path):
                raw_handle.write(line)

    join_genes_file = tools.fileOps.get_tmp_toil_file()
    join_genes_gp = tools.fileOps.get_tmp_toil_file()
    cmd = [['joingenes', '-f', raw_gtf_fofn, '-o', '/dev/stdout'],
           ['grep', '-P', '\tAUGUSTUS\t(exon|CDS|start_codon|stop_codon|tts|tss)\t'],
           ['sed', ' s/jg/augPB-/g']]
    tools.procOps.run_proc(cmd, stdout=join_genes_file)

    # passing the joingenes output through gtfToGenePred then genePredToGtf fixes the sort order for homGeneMapping
    cmd = ['gtfToGenePred', '-genePredExt', join_genes_file, join_genes_gp]
    tools.procOps.run_proc(cmd)
    cmd = ['genePredToGtf', 'file', join_genes_gp, '-utr', '-honorCdsStat', '-source=augustusPB', join_genes_file]
    tools.procOps.run_proc(cmd)

    joined_gtf_file_id = job.fileStore.writeGlobalFile(join_genes_file)
    raw_gtf_file_id = job.fileStore.writeGlobalFile(raw_gtf_file)
    joined_gp_file_id = job.fileStore.writeGlobalFile(join_genes_gp)
    return raw_gtf_file_id, joined_gtf_file_id, joined_gp_file_id
