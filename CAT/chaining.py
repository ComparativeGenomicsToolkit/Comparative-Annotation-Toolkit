""""
Toil program to generate UCSC chains and nets between two genomes in a HAL file.
"""
import os
import logging
import tools.hal
import tools.procOps
import tools.fileOps
import tools.toilInterface
from toil.job import Job


def chaining(args, toil_options):
    """entry point to this program"""
    j = Job.wrapJobFn(setup, args, toil_options)
    i = Job.Runner.startToil(j, toil_options)


def setup(job, args, toil_options):
    """
    Entry function for chaining cactus alignments
    :param args: argument dictionary
    :return:
    """
    input_file_ids = upload_files(job, args)
    chrom_sizes_local_path = job.fileStore.readGlobalFile(input_file_ids['sizes'])
    return_file_ids = []
    for i, l in enumerate(open(chrom_sizes_local_path)):
        chrom, size = l.split()
        j = job.addChildJobFn(chain_by_chromosome, args, chrom, size, input_file_ids, memory='8G')
        return_file_ids.append(j.rv())
    job.addFollowOnJobFn(merge, return_file_ids, args, toil_options)


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
    hal, target_two_bit, query_two_bit = download_files(job, input_file_ids, args)
    cmd = [["halLiftover", "--outPSL", hal, args['ref_genome'], bed_path, args['genome'], "/dev/stdout"],
           ["pslPosTarget", "/dev/stdin", "/dev/stdout"],
           ["axtChain", "-psl", "-verbose=0", "-linearGap=medium", "/dev/stdin", target_two_bit, query_two_bit, chain]]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(chain)


def merge(job, chain_files, args, toil_options):
    """
    Merge together chain files.
    :param chain_files: list of fileStore file_ids
    :param args: argument dictionary
    :param toil_options: argument namespace for Toil launching
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
    tools.toilInterface.export_to_master(tmp_chain_file, args['chain_file'], toil_options.master_ip)


def upload_files(job, args):
    """load files to jobStore"""
    hal_file_id = job.fileStore.writeGlobalFile(args['hal'])
    chrom_sizes_file_id = job.fileStore.writeGlobalFile(args['query_sizes'])
    query_two_bit_file_id = job.fileStore.writeGlobalFile(args['query_two_bit'])
    target_two_bit_file_id = job.fileStore.writeGlobalFile(args['target_two_bit'])
    return {'hal': hal_file_id, 'sizes': chrom_sizes_file_id, 'query_two_bit': query_two_bit_file_id,
            'target_two_bit': target_two_bit_file_id}


def download_files(job, input_file_ids, args):
    """download all files from fileStore for chaining"""
    work_dir = job.fileStore.getLocalTempDir()
    hal_local_path = os.path.join(work_dir, os.path.basename(args['hal']))
    job.fileStore.readGlobalFile(input_file_ids['hal'], hal_local_path)
    target_two_bit_local_path = os.path.join(work_dir, os.path.basename(args['target_two_bit']))
    job.fileStore.readGlobalFile(input_file_ids['target_two_bit'], target_two_bit_local_path)
    query_two_bit_local_path = os.path.join(work_dir, os.path.basename(args['query_two_bit']))
    job.fileStore.readGlobalFile(input_file_ids['query_two_bit'], query_two_bit_local_path)
    return hal_local_path, target_two_bit_local_path, query_two_bit_local_path
