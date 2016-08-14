"""
 file:    augustus_cgp.py
 descr.:  runs AugustusCGP on input HAL file
          optionally, a database with extrinsic evidence can be provided
          for parallel computing the HAL alignment is splitted into overlapping alignment chunks.
          Gene prediction chunks are merged with the auxiliary tool 'joingenes' from the
          Augustus package. The output is one gff file for each species in the clade
          (or the set of target genomes)

 authors: Stefanie Koenig, Ian Fiddes
 
  date    |  author         |  changes
 ---------|-----------------|------------------------------------------
 10.08.16 | Stefanie Koenig | creation of the file
"""

from toil.job import Job
from toil.common import Toil

import logging
import os

import tools.fileOps
import tools.procOps

def joinGenes(job, gffChunks):
    """
    uses the auxiliary tool 'joingenes' from the
    Augustus package to intelligently merge gene sets
    - removes duplicated Txs or truncated Txs that are contained in other Txs (trivial)
    - fixes truncated Txs at alignment boundaries,
      e.g. by merging them with other Txs (non trivial, introduces new Txs)
    Note: the 'grep' command removes transcript/gene features that are not
    supported by the GTF format. This is only a temporary fix and can be removed,
    once the gff files are processed downstream 'pythonicly'.
    """
    jg = tools.fileOps.get_tmp_file(tmp_dir=job.fileStore.getLocalTempDir())
    chunks = ",".join(map(job.fileStore.readGlobalFile,gffChunks))
    cmd = [['joingenes','-g', chunks, '-o', '/dev/stdout'],
           ['grep','-P','\tAUGUSTUS\t(exon|CDS|start_codon|stop_codon|tts|tss)\t']]
    tools.procOps.run_proc(cmd, stdout=jg)
    return job.fileStore.writeGlobalFile(jg)

def mergeGffs(job, args, gffChunks):
    """
    merges for each genome the gffChunks into one gff file
    """
    mergedGffs = {}
    for genome in args['fasta_files']:
        job.fileStore.logToMaster('Merging GFFs for {}'.format(genome), level=logging.INFO)
        # merge all gffChunks of one genome
        genome_gffChunks = [d[genome] for d in gffChunks]
        mergedGffs[genome] = job.addChildJobFn(joinGenes, genome_gffChunks).rv()
    return mergedGffs

def writeTree(job,input_file_ids):
    """
    writes a file with the phylogenetic tree in NEWICK format
    """
    hal = job.fileStore.readGlobalFile(input_file_ids['hal']) 
    cmd = ['halStats', '--tree', hal]
    tree = tools.fileOps.get_tmp_file(tmp_dir=job.fileStore.getLocalTempDir())
    tools.procOps.run_proc(cmd, stdout=tree)
    return job.fileStore.writeGlobalFile(tree)

def writeGenomeDirs(job,args):
    """
    writes a file with the location of the fasta files, e.g.

    galGal4 /path/to/genome/galGal4.fa
    hg38    /path/to/genome/hg38.fa
    mm10    /path/to/genome/mm10.fa
    rn6     /path/to/genome/rn6.fa
    ...
    """
    genomeDirs = tools.fileOps.get_tmp_file(tmp_dir=job.fileStore.getLocalTempDir())
    with open(genomeDirs, 'w') as outf:
        for genome in args['fasta_files']:
            outf.write(genome + '\t' + args['fasta_files'][genome] + '\n')
    return job.fileStore.writeGlobalFile(genomeDirs)

        
def hal2maf(job, input_file_ids, refGenome, chrom, start, chunkSize, genomic_region):
    """
    exports hal to maf on a genomic region specified by (genome, seq, start, len)
    """
    job.fileStore.logToMaster('Running hal2maf on {}'.format(genomic_region), level=logging.INFO)
    hal = job.fileStore.readGlobalFile(input_file_ids['hal']) 
    mafChunk = tools.fileOps.get_tmp_file(tmp_dir=job.fileStore.getLocalTempDir())
    cmd = ['hal2maf','--noAncestors', '--noDupes', '--refGenome', refGenome,
           '--refSequence', chrom, '--start', start, '--length', chunkSize, hal, mafChunk]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(mafChunk)

def cgp(job, tree, mafChunk, args, genomeDirs, genomic_region):
    """
    core function that runs AugustusCGP on one alignment chunk
    """
    job.fileStore.logToMaster('Running AugustusCGP on {}'.format(genomic_region), level=logging.INFO)

    work_dir = job.fileStore.getLocalTempDir()
    cgp_dir = os.path.join(work_dir, 'cgp') # where all the AugustusCGP output goes

    opt_param = []  # optional AugustusCGP parameters
    cgp_cfg = None
    if args['cgp_cfg'] is not None:
        # extrinsic config file
        # for now, let's assume RNA-Seq evidence whenever an extrinsic
        # config file is given, e.g. by default we turn on
        # --UTR=1
        # --allow_hinted_splicesites=atac
        cgp_cfg = job.fileStore.readGlobalFile(args['cgp_cfg'])
        opt_param += ['--extrinsicCfgFile={}'.format(cgp_cfg), 
                      '--dbhints=1', '--UTR=1', '--allow_hinted_splicesites=atac']
    if args['cgp_param'] is not None:
        # clade parameters from logistic regression
        # if not given the default parameters in the file
        # augustus/trunks/config/cgp/log_reg_parameters_default.cfg is used
        cgp_param = job.fileStore.readGlobalFile(args['cgp_param'])
        opt_param += ['--optCfgFile={}'.format(cgp_param)]

    cmd = ['augustus'] + opt_param + ['--species={}'.format(args['species']),
           '--treefile={}'.format(job.fileStore.readGlobalFile(tree)),
           '--alnfile={}'.format(job.fileStore.readGlobalFile(mafChunk)),
           '--dbaccess={}'.format(job.fileStore.readGlobalFile(args['hints_db'])),
           '--speciesfilenames={}'.format(job.fileStore.readGlobalFile(genomeDirs)),
           '--softmasking=1',
           '--exoncands=0',
           '--alternatives-from-evidence=0',
           '--/CompPred/logreg=on',
           '--printOEs=false',
           '--/CompPred/outdir={}'.format(cgp_dir)]

    aug_out = tools.procOps.call_proc_lines(cmd)
    return {genome: job.fileStore.writeGlobalFile(os.path.join(cgp_dir, genome + '.cgp.gff')) for genome in args['fasta_files']}

def setup(job, args, input_file_ids):
    """
    Entry function for running AugustusCGP.
    HAL alignment is converted to MAF format and splitted into overlapping
    alignment chunks for parallel computing. Each alignment chunk is one child process.
    Gene predictions on alignment chunks are subsequently merged into one gff for each species. 
    For merging of the gene sets, the auxiliary tool 'joingenes' from the Augustus tool package is used.
    """
    job.fileStore.logToMaster('Preparing input files for AugustusCGP', level=logging.INFO)
    
    # create a file with the locations of genome fasta file
    genomeDirs = writeGenomeDirs(job,args)
    
    # create a file with the phylogenetic tree in NEWICK format
    tree = writeTree(job,input_file_ids)

    # list of dicts, each storing all gffs for one alignment chunk
    # key: genome, value: file handle to gff
    gffChunks = []

    # calculate alignment chunks
    # TODO: do not split within genic regions of the reference genome
    
    # overlap length between two consecutive alignment chunks
    overlap = args['overlap']
    # length of alignment chunk with respect to the reference genome
    chunkSize = args['chunksize']

    aliChunks = [] # stores all alignment chunks as tuples [chrom, start, chunkSize]
    chromSizes = job.fileStore.readGlobalFile(input_file_ids['sizes'])
    for i, l in enumerate(open(chromSizes)):
        chrom, chromSize = l.split()
        start = 0
        while start + chunkSize < int(chromSize):
            aliChunks.append([chrom,start,chunkSize])
            start = start + chunkSize - overlap                # start of next alignment chunk
        aliChunks.append([chrom,start,int(chromSize)-start])   # last alignment chunk
    
    for a in aliChunks:
        genomic_region = '{}.{}:{}-{}'.format(args['ref_genome'],a[0],a[1],a[1]+a[2]-1) # string "genome.chrom:start-end"
        # export alignment chunks from hal to maf
        j = job.addChildJobFn(hal2maf, input_file_ids, args['ref_genome'], a[0], a[1], a[2], genomic_region)
        mafChunk = j.rv()
        # run AugustusCGP on alignment chunk
        gffChunks.append(j.addFollowOnJobFn(cgp,tree, mafChunk, args, genomeDirs, genomic_region).rv())
    
    # merge all gff files for alignment chunks to one gff for each species
    mergedGffs = job.addFollowOnJobFn(mergeGffs,args,gffChunks).rv()

    # TODO: assignment to parent genes

    return mergedGffs

def augustus_cgp(args, toil_options):
    """
    Main entry function for AugustusCGP toil pipeline
    :param args: dictionary of arguments from CAT
    :param toil_options: toil options Namespace object
    :return: dict with one gff file per species (keys=species, values=gff_file_handle)
    """
    with Toil(toil_options) as toil:
        if not toil.options.restart:
            fasta_files_id = toil.importFile('file:///' + args['hal'])
            hal_file_id = toil.importFile('file:///' + args['hal'])
            chrom_sizes_file_id = toil.importFile('file:///' + args['query_sizes'])
            input_file_ids = {'hal': hal_file_id, 'sizes': chrom_sizes_file_id}
            hints_db_file_id = toil.importFile('file:///' + args['hints_db'])
            if args['cgp_cfg'] is not None:
                cgp_cfg = toil.importFile('file:///' + args['cgp_cfg'])
            if args['cgp_param'] is not None:
                cgp_param = toil.importFile('file:///' + args['cgp_param'])

            job = Job.wrapJobFn(setup, args, input_file_ids)
            results = toil.start(job)
        else:
            results = toil.restart()
        for genome in results:
            tools.fileOps.ensure_file_dir(args['augustus_cgp_gtf'][genome])
            toil.exportFile(results[genome], 'file:///' + args['augustus_cgp_gtf'][genome])
