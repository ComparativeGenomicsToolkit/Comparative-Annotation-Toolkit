"""
Runs AugustusPB on a target genome
"""

import argparse
from toil.common import Toil
from toil.job import Job

import tools.bio
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
            gtf_file_id, df = toil.start(job)
        else:
            gtf_file_id, df = toil.restart()
        tools.fileOps.ensure_file_dir(args.augustus_pb_gtf)
        toil.exportFile(gtf_file_id, 'file://' + args.augustus_pb_gtf)
        return df


def setup(job, args, input_file_ids):
    """
    Entry function for running AugustusPB.
    The genome is chunked up and the resulting gene sets merged using joingenes.
    """
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.genome_fasta,
                                                                 prefix='genome', upper=False)

    predictions = []
    for chrom in genome_fasta:
        for start in xrange(0, len(genome_fasta[chrom]), args.pb_genome_chunksize - args.pb_genome_overlap):
            j = job.addChildJobFn(augustus_pb_chunk, args, input_file_ids, chrom, start, start + args.pb_genome_overlap)
            predictions.append(j.rv())

    # results contains a pair of [gff_file_id, dataframe] where the dataframe contains the alternative parental txs
    results = job.addFollowOnJobFn(join_genes, input_file_ids, predictions, memory='8G').rv()
    return results


def augustus_pb_chunk(job, args, input_file_ids, chrom, start, stop):
    """
    core function that runs AugustusPB on one genome chunk
    """
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids.genome_fasta,
                                                                 prefix='genome', upper=False)
    hints = job.fileStore.readGlobalFile(input_file_ids.hints_gff)
    pb_cfg = job.fileStore.readGlobalFile(input_file_ids.pb_cfg)
    tmp_fasta = tools.fileOps.get_tmp_toil_file()
    tools.bio.write_fasta(tmp_fasta, chrom, genome_fasta[chrom][start:stop])
    results = tools.fileOps.get_tmp_toil_file()

    cmd = ['augustus', '--UTR=1', '--softmasking=1', '--allow_hinted_splicesites=atac',
           '--alternatives-from-evidence=1',
           '--hintsfile={}'.format(hints),
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
    fofn = tools.fileOps.get_tmp_toil_file()
    with open(fofn, 'w') as outf:
        for chunk in predictions:
            local_path = job.fileStore.readGlobalFile(chunk)
            outf.write(local_path + '\n')

    jg = tools.fileOps.get_tmp_toil_file()
    cmd = [['joingenes', '-f', fofn, '-o', '/dev/stdout', '--alternatives'],
           ['grep', '-P', '\tAUGUSTUS\t(exon|CDS|start_codon|stop_codon|tts|tss)\t'],
           ['sed', ' s/jg/augPB-/g']]
    tools.procOps.run_proc(cmd, stdout=jg)
    joined_file_id = job.fileStore.writeGlobalFile(jg)
    j = job.addFollowOnJobFn(tools.parentGeneAssignment.assign_parents, input_file_ids.ref_db_path,
                             input_file_ids.filtered_tm_gp, input_file_ids.unfiltered_tm_gp,
                             joined_file_id, memory='8G')
    return j.rv()
