""""
Toil program to generate UCSC chains and nets between two genomes in a HAL file.
"""
import os
import logging
import tools.hal
import tools.procOps
import tools.fileOps
from toil.job import Job
from toil.common import Toil


def chaining(args, toil_options):
    """entry point to this program"""
    with Toil(toil_options) as toil:
        if not toil.options.restart:
            hal_file_id = toil.importFile('file:///' + args['hal'])
            chrom_sizes_file_id = toil.importFile('file:///' + args['query_sizes'])
            query_two_bit_file_id = toil.importFile('file:///' + args['query_two_bit'])
            target_two_bit_file_id = toil.importFile('file:///' + args['target_two_bit'])
            input_file_ids = {'hal': hal_file_id, 'sizes': chrom_sizes_file_id, 'query_two_bit': query_two_bit_file_id,
                              'target_two_bit': target_two_bit_file_id}
            job = Job.wrapJobFn(setup, args, input_file_ids)
            chain_file_id = toil.start(job)
        else:
            chain_file_id = toil.restart()
        tools.fileOps.ensure_file_dir(args['chain_file'])
        toil.exportFile(chain_file_id, 'file:///' + args['chain_file'])


def setup(job, args, input_file_ids):
    """
    Entry function for chaining cactus alignments
    :param args: argument dictionary
    :param input_file_ids: file ID dictionary of imported files
    :return: fileStore ID for output chain file
    """
    chrom_sizes = job.fileStore.readGlobalFile(input_file_ids['sizes'])
    return_file_ids = []
    for i, l in enumerate(open(chrom_sizes)):
        chrom, size = l.split()
        j = job.addChildJobFn(chain_by_chromosome, args, chrom, size, input_file_ids, memory='8G')
        return_file_ids.append(j.rv())
    return job.addFollowOnJobFn(merge, return_file_ids, args).rv()


def chain_by_chromosome(job, args, chrom, size, input_file_ids):
    """
    Chain alignments per-chromosome.
    :param args: argument dictionary
    :param chrom: chromosome name
    :param size: chromosome size
    :param input_file_ids: dict of file IDs in fileStore
    :return: chain file for this chromosome
    """
    job.fileStore.logToMaster('Beginning to chain chromosome {}-{}'.format(args['genome'], chrom),
                              level=logging.INFO)
    bed_path = tools.fileOps.get_tmp_file(tmp_dir=job.fileStore.getLocalTempDir())
    with open(bed_path, 'w') as outf:
        tools.fileOps.print_row(outf, [chrom, 0, size])
    chain = tools.fileOps.get_tmp_file(tmp_dir=job.fileStore.getLocalTempDir())
    # load files from jobStore
    work_dir = job.fileStore.getLocalTempDir()
    hal = os.path.join(work_dir, os.path.basename(args['hal']))
    job.fileStore.readGlobalFile(input_file_ids['hal'], hal)
    target_two_bit = os.path.join(work_dir, os.path.basename(args['target_two_bit']))
    job.fileStore.readGlobalFile(input_file_ids['target_two_bit'], target_two_bit)
    query_two_bit = os.path.join(work_dir, os.path.basename(args['query_two_bit']))
    job.fileStore.readGlobalFile(input_file_ids['query_two_bit'], query_two_bit)
    # execute liftover
    cmd = [["halLiftover", "--outPSL", hal, args['ref_genome'], bed_path, args['genome'], "/dev/stdout"],
           ["pslPosTarget", "/dev/stdin", "/dev/stdout"],
           ["axtChain", "-psl", "-verbose=0", "-linearGap=medium", "/dev/stdin", target_two_bit, query_two_bit, chain]]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(chain)


def merge(job, chain_files, args):
    """
    Merge together chain files.
    :param chain_files: list of fileStore file_ids
    :param args: argument dictionary
    :return:
    """
    job.fileStore.logToMaster('Merging chains for {}'.format(args['genome']), level=logging.INFO)
    fofn = tools.fileOps.get_tmp_file(tmp_dir=job.fileStore.getLocalTempDir())
    with open(fofn, 'w') as outf:
        for file_id in chain_files:
            local_path = job.fileStore.readGlobalFile(file_id)
            outf.write(local_path + '\n')
    cmd = ['chainMergeSort', '-inputList={}'.format(fofn), '-tempDir={}/'.format(job.fileStore.getLocalTempDir())]
    tmp_chain_file = tools.fileOps.get_tmp_file(tmp_dir=job.fileStore.getLocalTempDir())
    tools.procOps.run_proc(cmd, stdout=tmp_chain_file)
    tmp_chain_file_id = job.fileStore.writeGlobalFile(tmp_chain_file)
    return tmp_chain_file_id
