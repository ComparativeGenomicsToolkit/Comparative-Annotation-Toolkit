""""
Toil program to generate UCSC chains and nets between two genomes in a HAL file.
"""
import logging
import tools.hal
import tools.procOps
import tools.fileOps
from toil.job import Job


def chain_by_chromosome(job, args, hal, chrom, size):
    """
    Chain alignments per-chromosome.
    :param job: job
    :param hal: hal alignment file
    :param args: argument dictionary
    :param chrom: chromosome name
    :param size: chromosome size
    :return: chain file for this chromosome
    """
    job.fileStore.logToMaster('Beginning to chain chromosome {}-{}'.format(args['genome'], chrom),
                              level=logging.INFO)
    bed_path = tools.fileOps.get_tmp_file(tmp_dir=job.fileStore.getLocalTempDir())
    with open(bed_path, 'w') as outf:
        tools.fileOps.print_row(outf, chrom, 0, size)
    chain = tools.fileOps.get_tmp_file(tmp_dir=job.fileStore.getLocalTempDir())
    cmd = [["halLiftover", "--outPSL", hal, args['ref_genome'], bed_path, args['genome'], "/dev/stdout"],
           ["pslPosTarget", "/dev/stdin", "/dev/stdout"],
           ["axtChain", "-psl", "-verbose=0", "-linearGap=medium", "/dev/stdin", args['target_two_bit'],
            args['query_two_bit'], chain]]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(chain)


def merge(job, chain_files, args):
    """
    Merge together chain files.
    :param job: job
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
    tools.fileOps.ensure_file_dir(args['chain_file'])
    cmd = ['chainMergeSort', '-inputList={}'.format(fofn), '-tempDir={}/'.format(job.fileStore.getLocalTempDir())]
    tools.procOps.run_proc(cmd, stdout=args['chain_file'])


def setup(job, args, hal):
    """
    Entry function for chaining cactus alignments
    :param hal: HAL file
    :param args: argument dictionary
    :return:
    """
    file_ids = []
    for i, l in enumerate(open(args['query_sizes'])):
        chrom, size = l.split()
        j = job.addChildJobFn(chain_by_chromosome, args, hal, chrom, size, memory='8G')
        file_ids.append(j.rv())
    job.addFollowOnJobFn(merge, file_ids, args)


def chaining(args, hal, toil_options):
    j = Job.wrapJobFn(setup, args, hal)
    i = Job.Runner.startToil(j, toil_options)
