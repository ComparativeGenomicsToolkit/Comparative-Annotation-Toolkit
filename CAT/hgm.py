"""
 file:    hgm.py
 descr.:  runs homGeneMapping from the Augustus package
          for cross-species comparison of gene sets. HomGeneMapping uses
          the intron hints from the database to retrieve on a per-transcript
          basis information, which introns have RNA-Seq splice junctions (SJ) supported in which of the
          input genomes. It essentially adds the "hgm_info" string to the last column of the gtf, 

          hgm_info "16E-27,0E-13,1,2E-13,3E-1,4,5,6E-48,7E-30,8E-1,9,10E-1,11,12E-19,13E-27,14E-46,15E-6,17E-4";
          
          that encodes genome name, type of evidence and multiplicity, e.g.
          in the example above, the intron has RNA-Seq SJ support in species 16 (with mult=26),
          in species 0 (with mult=13), in species 2 (with mult=13), in species 3 (with mult=1), etc.
          The header in the gtf, gives a list of species numbers and corresponding names, e.g.
          
          # 0     129S1_SvImJ
          # 1     AKR_J
          # 2     A_J
          ...
          # 17    WSB_EiJ

 authors: Stefanie Koenig, Ian Fiddes
 
  date    |  author         |  changes
 ---------|-----------------|------------------------------------------
 31.08.16 | Stefanie Koenig | creation of the file
"""

from toil.job import Job
from toil.common import Toil

import logging
import os
import collections
import multiprocessing

import tools.fileOps
import tools.procOps

###
# hgm pipeline section
###


def hgm(args, toil_options):
    """
    Main entry function for hgm toil pipeline
    :param args: dictionary of arguments from CAT
    :param toil_options: toil options Namespace object
    :return: a dictionary with one gtf file per genome
    """
    # we decide to either use toil_options.maxCores or the # of CPU on the leader node
    # this could lock up if the leader has more power than the children
    # we pass this along in the args for homGeneMapping to make use of
    num_cpu = min(multiprocessing.cpu_count(), toil_options.maxCores)
    with Toil(toil_options) as toil:
        if not toil.options.restart:
            hal_file_id = toil.importFile('file://' + args['hal'])
            hints_db_file_id = toil.importFile('file://' + args['hints_db'])
            gtf_file_ids = {genome: toil.importFile('file://' + gtf) for genome, gtf in args['in_gtf'].iteritems()}
            input_file_ids = {'hal': hal_file_id, 'hints_db': hints_db_file_id, 'gtfs': gtf_file_ids}
            job = Job.wrapJobFn(setup, args, input_file_ids, num_cpu, cores=num_cpu)
            results = toil.start(job)
        else:
            results = toil.restart()
        for genome in results:
            tools.fileOps.ensure_file_dir(args['hgm_gtf'][genome])
            toil.exportFile(results[genome], 'file://' + args['hgm_gtf'][genome])


def writeGtfFofn(job, gtf_file_ids):
    """
    writes a file with the location of the input gtf files, e.g.
    galGal4 /path/to/gene_set/galGal4.gtf
    hg38    /path/to/gene_set/hg38.gtf
    mm10    /path/to/gene_set/mm10.gtf
    rn6     /path/to/gene_set/rn6.gtf
    ...
    These files are loaded from the fileStore
    """
    gtfFofn = tools.fileOps.get_tmp_toil_file()
    with open(gtfFofn, 'w') as outf:
        for genome, file_id in gtf_file_ids.iteritems():
            local_path = job.fileStore.readGlobalFile(file_id)
            tools.fileOps.print_row(outf, [genome, local_path])
    return gtfFofn


def setup(job, args, input_file_ids, num_cpu):
    """
    runs homGeneMapping on a set of input gtf files from different genomes and
    returns one gtf file for each genome
    containing the "hgm_info" string in the last column
    """
    job.fileStore.logToMaster('Running homGeneMapping with {} cores'.format(num_cpu), level=logging.INFO)
    gtfFofn = writeGtfFofn(job, input_file_ids['gtfs'])

    cmd = ['homGeneMapping',
           '--halfile={}'.format(job.fileStore.readGlobalFile(input_file_ids['hal'])),
           '--dbaccess={}'.format(job.fileStore.readGlobalFile(input_file_ids['hints_db'])),
           '--gtfs={}'.format(gtfFofn),
           '--cpus={}'.format(num_cpu),
           '--outdir={}'.format(os.getcwd())]
    tools.procOps.run_proc(cmd)
    return {genome: job.fileStore.writeGlobalFile(genome + '.gtf') for genome in args['genomes']}    


def parse_hgm_gtf(hgm_out):
    """
    parses the hgm output gtfs and creates for each transcript a string with the intron support counts
    For now, we just count for each of the introns in the transcript, the number of species
    in which it has RNA-Seq SJ support (number of "E" in the 'hgm_info' string).
    But this can be changed later on, e.g. using also the multiplicities, or the presence of the introng
    in (one of) the reference annotation(s) (if "M" is in the 'hgm_info' string, then it is an annotated intron)
    """
    d = collections.defaultdict(list)
    
    with open(hgm_out, 'r') as infile:
        # get last column of all intron lines
        intron_lines = [i.strip().split('\t')[-1] for i in infile if "\tintron\t" in i]
        for attributes in intron_lines:
            txid = None
            count = 0
            for a in attributes.split(';'):
                if "hgm_info" in a:                
                    count = a.count("E")  # count number of occurrences of 'E'
                if "transcript_id" in a:
                    txid = a.split()[-1].strip('"')  # parse transcript id
            if txid is None:
                raise RuntimeError("Internal error in parse_hgm_gtf. Missing transcript_id in file {}".format(hgm_out))
            d[txid].append(count)
    return {k: ','.join(map(str, v)) for k, v in d.items()}  # convert list of intron counts to comma-separated string
