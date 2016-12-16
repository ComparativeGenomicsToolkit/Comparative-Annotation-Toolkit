"""
 file:    augustus_cgp.py
 descr.:  runs AugustusCGP on input HAL file
          optionally, a database with extrinsic evidence can be provided
          for parallel computing the HAL alignment is split into overlapping alignment chunks.
          Gene prediction chunks are merged with the auxiliary tool 'joingenes' from the
          Augustus package. The output is one gff file for each species in the clade
          (or the set of target genomes)

 authors: Stefanie Koenig, Ian Fiddes
 
  date    |  author         |  changes
 ---------|-----------------|------------------------------------------
 10.08.16 | Stefanie Koenig | creation of the file
"""

import argparse
import os

from toil.common import Toil
from toil.job import Job

import tools.dataOps
import tools.fileOps
import tools.intervals
import tools.procOps
import tools.sqlInterface
import tools.transcripts
import tools.parentGeneAssignment


def augustus_cgp(args, toil_options):
    """
    Main entry function for AugustusCGP toil pipeline
    :param args: dictionary of arguments from CAT
    :param toil_options: toil options Namespace object
    :return:
    """
    with Toil(toil_options) as toil:
        if not toil.options.restart:
            input_file_ids = argparse.Namespace()
            input_file_ids.hal = toil.importFile('file://' + args.hal)
            input_file_ids.chrom_sizes = toil.importFile('file://' + args.query_sizes)
            input_file_ids.hints_db = toil.importFile('file://' + args.hints_db)
            input_file_ids.cgp_param = toil.importFile('file://' + args.cgp_param)
            input_file_ids.ref_db_path = toil.importFile('file://' + args.ref_db_path)
            input_file_ids.fasta = {genome: toil.importFile('file://' + fasta)
                                    for genome, fasta in args.fasta_files.iteritems()}
            input_file_ids.filtered_tm_gps = {genome: toil.importFile('file://' + tm_gp)
                                              for genome, tm_gp in args.filtered_tm_gps.iteritems()}
            input_file_ids.unfiltered_tm_gps = {genome: toil.importFile('file://' + tm_gp)
                                                for genome, tm_gp in args.unfiltered_tm_gps.iteritems()}
            input_file_ids.cgp_cfg = toil.importFile('file://' + args.cgp_cfg)
            job = Job.wrapJobFn(setup, args, input_file_ids, memory='8G')
            results = toil.start(job)
        else:
            results = toil.restart()
        dataframes = []
        for genome, (gtf_file_id, df) in results.iteritems():
            tools.fileOps.ensure_file_dir(args.augustus_cgp_gtf[genome])
            toil.exportFile(gtf_file_id, 'file://' + args.augustus_cgp_gtf[genome])
            dataframes.append([genome, df])
        return dataframes


def setup(job, args, input_file_ids):
    """
    Entry function for running AugustusCGP.
    HAL alignment is converted to MAF format and splitted into overlapping
    alignment chunks for parallel computing. Each alignment chunk is one child process.
    Gene predictions on alignment chunks are subsequently merged into one gff for each species.
    For merging of the gene sets, the auxiliary tool 'joingenes' from the Augustus tool package is used.
    """
    # create a file with the phylogenetic tree in NEWICK format
    tree = writeTree(job, input_file_ids)

    # list of dicts, each storing all gffs for one alignment chunk
    # key: genome, value: file handle to gff
    gffChunks = []

    # calculate alignment chunks
    # TODO: do not split within genic regions of the reference genome

    # overlap length between two consecutive alignment chunks
    overlap = args.overlap
    # length of alignment chunk with respect to the reference genome
    chunkSize = args.chunksize

    aliChunks = []  # stores all alignment chunks as tuples [chrom, start, chunkSize]
    chromSizes = job.fileStore.readGlobalFile(input_file_ids.chrom_sizes)
    for chrom, chromSize in tools.fileOps.iter_lines(chromSizes):
        start = 0
        while start + chunkSize < int(chromSize):
            aliChunks.append([chrom, start, chunkSize])
            start = start + chunkSize - overlap                  # start of next alignment chunk
        aliChunks.append([chrom, start, int(chromSize) - start])   # last alignment chunk

    for chrom, start, end in aliChunks:
        # export alignment chunks from hal to maf
        j = job.addChildJobFn(hal2maf, input_file_ids, args.ref_genome, chrom, start, end, memory='8G')
        mafChunk = j.rv()
        # run AugustusCGP on alignment chunk
        cgp_job = j.addFollowOnJobFn(cgp, tree, mafChunk, args, input_file_ids, memory='8G')
        gffChunk = cgp_job.rv()
        gffChunks.append(gffChunk)

    # merge all gff files for alignment chunks to one gff for each species
    # results contains pairs of [gff_file_id, dataframe] where the dataframe contains the alternative parental txs
    results = job.addFollowOnJobFn(merge_results, args, input_file_ids, gffChunks, memory='8G').rv()
    return results


def hal2maf(job, input_file_ids, refGenome, chrom, start, chunkSize):
    """
    exports hal to maf on a genomic region specified by (genome, seq, start, len)
    """
    hal = job.fileStore.readGlobalFile(input_file_ids.hal)
    mafChunk = tools.fileOps.get_tmp_toil_file()
    cmd = ['hal2maf', '--noAncestors', '--noDupes', '--refGenome', refGenome,
           '--refSequence', chrom, '--start', start, '--length', chunkSize, hal, mafChunk]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(mafChunk)


def cgp(job, tree, mafChunk, args, input_file_ids):
    """
    core function that runs AugustusCGP on one alignment chunk
    """
    genomeFofn = writeGenomeFofn(job, input_file_ids.fasta)
    cgp_cfg = job.fileStore.readGlobalFile(input_file_ids.cgp_cfg)

    cmd = ['augustus', '--dbhints=1', '--UTR=1', '--allow_hinted_splicesites=atac',
           '--extrinsicCfgFile={}'.format(cgp_cfg),
           '--species={}'.format(args.species),
           '--treefile={}'.format(job.fileStore.readGlobalFile(tree)),
           '--alnfile={}'.format(job.fileStore.readGlobalFile(mafChunk)),
           '--dbaccess={}'.format(job.fileStore.readGlobalFile(input_file_ids.hints_db)),
           '--speciesfilenames={}'.format(genomeFofn),
           '--softmasking=1',
           '--exoncands=0',
           '--alternatives-from-evidence=0',
           '--/CompPred/logreg=on',
           '--printOEs=false',
           '--/CompPred/outdir={}'.format(os.getcwd()),
           '--optCfgFile={}'.format(job.fileStore.readGlobalFile(input_file_ids.cgp_param))]
    tools.procOps.run_proc(cmd)
    return {genome: job.fileStore.writeGlobalFile(genome + '.cgp.gff') for genome in args.genomes}


def merge_results(job, args, input_file_ids, gffChunks):
    """
    Merges the results using joinGenes. The results have parental genes assigned.
    """
    results = {}
    for genome in args.genomes:
        # merge all gffChunks of one genome
        genome_gffChunks = [d[genome] for d in gffChunks]
        j = job.addChildJobFn(join_genes, genome, input_file_ids, genome_gffChunks, memory='8G')
        results[genome] = j.rv()
    return results


def join_genes(job, genome, input_file_ids, gffChunks):
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
        for chunk in gffChunks:
            local_path = job.fileStore.readGlobalFile(chunk)
            outf.write(local_path + '\n')

    jg = tools.fileOps.get_tmp_toil_file()
    cmd = [['joingenes', '-f', fofn, '-o', '/dev/stdout'],
           ['grep', '-P', '\tAUGUSTUS\t(exon|CDS|start_codon|stop_codon|tts|tss)\t'],
           ['sed', ' s/jg/augCGP_/g']]
    tools.procOps.run_proc(cmd, stdout=jg)
    joined_file_id = job.fileStore.writeGlobalFile(jg)
    j = job.addFollowOnJobFn(tools.parentGeneAssignment.assign_parents, input_file_ids.ref_db_path,
                             input_file_ids.filtered_tm_gps[genome], input_file_ids.unfiltered_tm_gps[genome],
                             joined_file_id, memory='8G')
    return j.rv()


###
# Accessory functions
###


def writeTree(job,input_file_ids):
    """
    writes a file with the phylogenetic tree in NEWICK format
    """
    hal = job.fileStore.readGlobalFile(input_file_ids.hal) 
    cmd = ['halStats', '--tree', hal]
    tree = tools.fileOps.get_tmp_toil_file()
    tools.procOps.run_proc(cmd, stdout=tree)
    return job.fileStore.writeGlobalFile(tree)


def writeGenomeFofn(job, fasta_file_ids):
    """
    writes a file with the location of the fasta files, e.g.

    galGal4 /path/to/genome/galGal4.fa
    hg38    /path/to/genome/hg38.fa
    mm10    /path/to/genome/mm10.fa
    rn6     /path/to/genome/rn6.fa
    ...

    These files are loaded from the fileStore
    """
    genomeFofn = tools.fileOps.get_tmp_toil_file()
    with open(genomeFofn, 'w') as outf:
        for genome, file_id in fasta_file_ids.iteritems():
            local_path = job.fileStore.readGlobalFile(file_id)
            tools.fileOps.print_row(outf, [genome, local_path])
    return genomeFofn


