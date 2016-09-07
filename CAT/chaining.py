""""
Toil program to generate UCSC chains and nets between two genomes in a HAL file.
"""
import logging
import collections
import tools.hal
import tools.procOps
import tools.fileOps
from toil.job import Job
from toil.common import Toil


def chaining(args, toil_options):
    """entry point to this program"""
    with Toil(toil_options) as toil:
        if not toil.options.restart:
            target_two_bit_file_ids = {genome: toil.importFile('file://' + f) for genome, f
                                       in args['target_two_bits'].iteritems()}
            input_file_ids = {'hal': toil.importFile('file://' + args['hal']),
                              'query_sizes': toil.importFile('file://' + args['query_sizes']),
                              'query_two_bit': toil.importFile('file://' + args['query_two_bit']),
                              'target_two_bits': target_two_bit_file_ids}
            job = Job.wrapJobFn(setup, args, input_file_ids)
            chain_file_ids = toil.start(job)
        else:
            chain_file_ids = toil.restart()
        for chain_file, chain_file_id in chain_file_ids.iteritems():
            tools.fileOps.ensure_file_dir(chain_file)
            toil.exportFile(chain_file_id, 'file://' + chain_file)


def setup(job, args, input_file_ids):
    """
    Entry function for chaining cactus alignments
    :param args: argument dictionary
    :param input_file_ids: file ID dictionary of imported files
    :return: fileStore ID for output chain file
    """
    chrom_sizes = job.fileStore.readGlobalFile(input_file_ids['query_sizes'])
    tmp_chain_file_ids = collections.defaultdict(list)
    for i, l in enumerate(open(chrom_sizes)):
        chrom, size = l.split()
        for target_genome, target_two_bit_file_id in input_file_ids['target_two_bits'].iteritems():
            j = job.addChildJobFn(chain_by_chromosome, args, chrom, size, input_file_ids, target_genome,
                                  target_two_bit_file_id, memory='8G')
            tmp_chain_file_ids[target_genome].append(j.rv())
    return_file_ids = {}
    for genome, chain_file in args['chain_files'].iteritems():
        chain_files = tmp_chain_file_ids[genome]
        j = job.addFollowOnJobFn(merge, chain_files, genome)
        return_file_ids[chain_file] = j.rv()
    return return_file_ids


def chain_by_chromosome(job, args, chrom, size, input_file_ids, target_genome, target_two_bit_file_id):
    """
    Chain alignments per-chromosome.
    :param args: argument dictionary
    :param chrom: chromosome name
    :param size: chromosome size
    :param input_file_ids: dict of file IDs in fileStore
    :param target_genome: the genome we are analyzing here
    :param target_two_bit_file_id: the file ID for the twobit file for target_genome
    :return: chain file for this chromosome
    """
    job.fileStore.logToMaster('Beginning to chain chromosome {}-{}'.format(target_genome, chrom),
                              level=logging.INFO)
    bed_path = tools.fileOps.get_tmp_toil_file()
    with open(bed_path, 'w') as outf:
        tools.fileOps.print_row(outf, [chrom, 0, size])
    chain = tools.fileOps.get_tmp_toil_file()
    # load files from jobStore
    hal = job.fileStore.readGlobalFile(input_file_ids['hal'])
    target_two_bit = job.fileStore.readGlobalFile(target_two_bit_file_id)
    query_two_bit = job.fileStore.readGlobalFile(input_file_ids['query_two_bit'])
    # execute liftover
    cmd = [['halLiftover', '--outPSL', hal, args['ref_genome'], bed_path, target_genome, '/dev/stdout'],
           ['pslPosTarget', '/dev/stdin', '/dev/stdout'],
           ['axtChain', '-psl', '-verbose=0', '-linearGap=medium', '/dev/stdin', target_two_bit, query_two_bit, chain]]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(chain)


def merge(job, chain_files, genome):
    """
    Merge together chain files.
    :param chain_files: list of fileStore file_ids
    :param genome: genome being combined
    :return:
    """
    job.fileStore.logToMaster('Merging chains for {}'.format(genome), level=logging.INFO)
    fofn = tools.fileOps.get_tmp_toil_file()
    with open(fofn, 'w') as outf:
        for file_id in chain_files:
            local_path = job.fileStore.readGlobalFile(file_id)
            outf.write(local_path + '\n')
    cmd = ['chainMergeSort', '-inputList={}'.format(fofn), '-tempDir={}/'.format(job.fileStore.getLocalTempDir())]
    tmp_chain_file = tools.fileOps.get_tmp_toil_file()
    tools.procOps.run_proc(cmd, stdout=tmp_chain_file)
    tmp_chain_file_id = job.fileStore.writeGlobalFile(tmp_chain_file)
    return tmp_chain_file_id
