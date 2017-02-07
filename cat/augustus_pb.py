"""
Runs AugustusPB on a target genome
"""
import argparse
import collections

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
            job = Job.wrapJobFn(setup, args, input_file_ids, memory='8G', disk='2G')
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
    sizes_file = job.fileStore.readGlobalFile(input_file_ids.chrom_sizes)
    disk_usage = tools.toilInterface.find_total_disk_usage([input_file_ids.genome_fasta, input_file_ids.hints_gff])

    predictions = []
    for chrom, size in tools.fileOps.iter_lines(sizes_file):
        j = job.addChildJobFn(augustus_pb_chunk, args, input_file_ids, chrom, memory='32G', disk=disk_usage)
        predictions.append(j.rv())

    # results contains a 3 member tuple of [raw_gtf_file_id, joined_gtf_file_id, joined_gp_file_id]
    results = job.addFollowOnJobFn(combine_results, predictions, memory='2G', disk='8G').rv()
    return results


def augustus_pb_chunk(job, args, input_file_ids, chrom):
    """
    core function that runs AugustusPB on one genome chunk
    """
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.genome_fasta,
                                                                 prefix='genome', upper=False)
    hints = job.fileStore.readGlobalFile(input_file_ids.hints_gff)

    # slice out only the relevant hints, upgrade PB priority
    hints_subset = tools.fileOps.get_tmp_toil_file()
    cmd = [['awk', '($1 == "{}") {{print $0}}'.format(chrom)],
           ['sed', 's/pri=4;src=PB/pri=5;src=PB/g']]
    tools.procOps.run_proc(cmd, stdin=hints, stdout=hints_subset)

    pb_cfg = job.fileStore.readGlobalFile(input_file_ids.pb_cfg)
    tmp_fasta = tools.fileOps.get_tmp_toil_file()
    tools.bio.write_fasta(tmp_fasta, chrom, str(genome_fasta[chrom]))
    results = tools.fileOps.get_tmp_toil_file()

    cmd = ['augustus', '--UTR=1', '--softmasking=1', '--allow_hinted_splicesites=atac',
           '--alternatives-from-evidence=1',
           '--hintsfile={}'.format(hints_subset),
           '--extrinsicCfgFile={}'.format(pb_cfg),
           '--species={}'.format(args.species),
           tmp_fasta]
    tools.procOps.run_proc(cmd, stdout=results)
    return job.fileStore.writeGlobalFile(results)


def combine_results(job, predictions):
    """
    Combines the results. We can't use joingenes because it can't handle the alternative transcripts.
    """
    # first write the full raw results
    raw_gtf_file = tools.fileOps.get_tmp_toil_file()
    with open(raw_gtf_file, 'w') as raw_handle:
        for chunk in predictions:
            local_path = job.fileStore.readGlobalFile(chunk)
            for line in open(local_path):
                raw_handle.write(line)

    # filter for relevant entries
    tmp_gtf_file = tools.fileOps.get_tmp_toil_file()
    cmd = ['grep', '-P', '\tAUGUSTUS\t(exon|CDS|start_codon|stop_codon|tts|tss)\t', raw_gtf_file]
    tools.procOps.run_proc(cmd, stdout=tmp_gtf_file)

    # use gtfToGenePred to convert this, and make it easier to convert names
    tmp_gp = tools.fileOps.get_tmp_toil_file()
    cmd = ['gtfToGenePred', '-genePredExt', tmp_gtf_file, tmp_gp]
    tools.procOps.run_proc(cmd)

    # create a dict mapping chrom -> gene_id -> tx objs
    chrom_gene_map = collections.defaultdict(lambda: collections.defaultdict(list))
    for tx in tools.transcripts.gene_pred_iterator(tmp_gp):
        chrom_gene_map[tx.chromosome][tx.name2].append(tx)

    # now fix the names
    gene_count = 0
    gps = []
    for chrom, gene_dict in chrom_gene_map.iteritems():
        for tx_list in gene_dict.itervalues():
            tx_count = 0
            gene_count += 1
            gene_id = 'augPB-g{}'.format(gene_count)
            for tx_obj in tx_list:
                tx_count += 1
                tx_id = 'augPB-g{}.t{}'.format(gene_count, tx_count)
                gps.append(tx_obj.get_gene_pred(name=tx_id, name2=gene_id))

    # write the results back as genePred
    joined_gp_file = tools.fileOps.get_tmp_toil_file()
    tools.fileOps.print_rows(joined_gp_file, gps)

    # convert back to GTF for final output
    join_genes_file = tools.fileOps.get_tmp_toil_file()
    cmd = ['genePredToGtf', 'file', joined_gp_file, '-utr', '-honorCdsStat', '-source=augustusPB', join_genes_file]
    tools.procOps.run_proc(cmd)

    joined_gtf_file_id = job.fileStore.writeGlobalFile(join_genes_file)
    raw_gtf_file_id = job.fileStore.writeGlobalFile(raw_gtf_file)
    joined_gp_file_id = job.fileStore.writeGlobalFile(joined_gp_file)
    return raw_gtf_file_id, joined_gtf_file_id, joined_gp_file_id
