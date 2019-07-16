"""
 file:    augustus_cgp.py
 descr.:  runs AugustusCGP on input HAL file
          optionally, a database with extrinsic evidence can be provided
          for parallel computing the HAL alignment is split into overlapping alignment chunks.
          Gene prediction chunks are merged with the auxiliary tool 'joingenes' from the
          Augustus package. The output is one gff file for each species in the clade
          (or the set of target genomes)

 authors: Stefanie Koenig, Ian Fiddes
"""

import argparse
import sqlalchemy
import os
import collections
import random

from toil.fileStore import FileID
from toil.common import Toil
from toil.job import Job

import tools.misc
import tools.toilInterface
import tools.dataOps
import tools.fileOps
import tools.intervals
import tools.procOps
import tools.sqlInterface
import tools.transcripts
import tools.hintsDatabaseInterface


def augustus_cgp(args, toil_options):
    """
    Main entry function for AugustusCGP toil pipeline
    :param args: dictionary of arguments from CAT
    :param toil_options: toil options Namespace object
    :return:
    """
    with Toil(toil_options) as t:
        if not t.options.restart:
            input_file_ids = argparse.Namespace()
            input_file_ids.hal = FileID.forPath(t.importFile('file://' + args.hal), args.hal)
            input_file_ids.chrom_sizes = FileID.forPath(t.importFile('file://' + args.query_sizes), args.query_sizes)
            input_file_ids.hints_db = FileID.forPath(t.importFile('file://' + args.hints_db), args.hints_db)
            if args.cgp_param is not None:
                input_file_ids.cgp_param = FileID.forPath(t.importFile('file://' + args.cgp_param), args.cgp_param)
            else:
                input_file_ids.cgp_param = None
                input_file_ids.gtf = FileID.forPath(t.importFile('file://' + args.gtf), args.gtf)
            input_file_ids.cgp_cfg = FileID.forPath(t.importFile('file://' + args.cgp_cfg), args.cgp_cfg)
            input_file_ids.fasta = {genome: FileID.forPath(t.importFile('file://' + fasta), fasta)
                                    for genome, fasta in args.fasta_files.iteritems()}
            du = tools.toilInterface.find_total_disk_usage([input_file_ids.hints_db], buffer='4G')
            job = Job.wrapJobFn(setup, args, input_file_ids, memory='8G', disk=du)
            results, stdout_file_ids, param_file_id = t.start(job)
        else:
            results, stdout_file_ids, param_file_id = t.restart()
        tools.fileOps.ensure_file_dir(args.stdout_file)
        with open(args.stdout_file, 'w') as outf, tools.fileOps.TemporaryFilePath() as tmp:
            for (chrom, start, chunksize), stdout_file in stdout_file_ids.iteritems():
                outf.write('## BEGIN CHUNK chrom: {} start: {} chunksize: {}\n'.format(chrom, start, chunksize))
                t.exportFile(stdout_file, 'file://' + tmp)
                for l in open(tmp):
                    outf.write(l)
        for genome, (raw_gtf_file_id, joined_gtf_file_id, joined_gp_file_id) in results.iteritems():
            tools.fileOps.ensure_file_dir(args.augustus_cgp_raw_gtf[genome])
            t.exportFile(raw_gtf_file_id, 'file://' + args.augustus_cgp_raw_gtf[genome])
            t.exportFile(joined_gtf_file_id, 'file://' + args.augustus_cgp_gtf[genome])
            t.exportFile(joined_gp_file_id, 'file://' + args.augustus_cgp_gp[genome])
        if args.cgp_param is None:
            t.exportFile(param_file_id, 'file://' + args.param_out_path)


def setup(job, args, input_file_ids):
    """
    Entry function for running AugustusCGP. Will first determine if we need to perform training, and do so.
    HAL alignment is converted to MAF format and split into overlapping
    alignment chunks for parallel computing. Each alignment chunk is one child process.
    Gene predictions on alignment chunks are subsequently merged into one gff for each species.
    For merging of the gene sets, the auxiliary tool 'joingenes' from the Augustus tool package is used.
    """
    # create a file with the phylogenetic tree in NEWICK format
    tree = write_tree(job, input_file_ids)

    # construct all MAF chunks
    chrom_sizes = job.fileStore.readGlobalFile(input_file_ids.chrom_sizes)
    # 4G buffer for MAF chunk, should be more than enough (famous last words)
    hal2maf_usage = tools.toilInterface.find_total_disk_usage(input_file_ids.hal)

    # TODO: do not split within genic regions of the reference genome
    maf_chunks = []  # list of lists [chrom, start, chunksize, fileID]
    for chrom, chrom_size in tools.fileOps.iter_lines(chrom_sizes):
        chrom_size = int(chrom_size)
        for start in xrange(0, chrom_size, args.chunksize - args.overlap):
            chunksize = args.chunksize if start + args.chunksize <= chrom_size else chrom_size - start
            j = job.addChildJobFn(hal2maf, input_file_ids, args.genomes, args.ref_genome, args.annotate_ancestors,
                                  chrom, start, chunksize, memory='8G', disk=hal2maf_usage)
            maf_chunks.append([chrom, start, chunksize, j.rv()])

    # if we have no params, time to train
    if input_file_ids.cgp_param is None:
        du = tools.toilInterface.find_total_disk_usage([input_file_ids.hints_db], buffer='40G')
        results = job.addFollowOnJobFn(cgp_training_wrapper, maf_chunks, tree, args, input_file_ids, memory='8G',
                                       disk=du).rv()
    else:
        results = job.addFollowOnJobFn(cgp_wrapper, maf_chunks, tree, args, input_file_ids, disk='4G').rv()
    return results


def cgp_training_wrapper(job, maf_chunks, tree, args, input_file_ids):
    """
    Wrapper around calling augustusCGP for training on a subset of the alignment.
    :param maf_chunks: List of lists [chrom, start, chunksize, fileID]
    :param tree: fileID for tree from write_tree()
    :param args: input arguments
    :param input_file_ids: input file IDs
    :return: output of merge_results()
    """
    hints_db = job.fileStore.readGlobalFile(input_file_ids.hints_db)

    # load hints database and seqnr information
    speciesnames, seqnames, hints, featuretypes, session = tools.hintsDatabaseInterface.reflect_hints_db(hints_db)
    speciesid = session.query(speciesnames.speciesid).filter_by(speciesname=args.ref_genome)
    seqs = {x.seqname: x.seqnr for x in session.query(seqnames).filter_by(speciesid=speciesid)}

    # begin selecting random intervals
    selected_intervals = []
    seen_exons = 0
    for chrom, start, chunksize, maf_chunk in random.sample(maf_chunks, len(maf_chunks)):
        seqnr = seqs[chrom]
        query = session.query(hints).filter(
                sqlalchemy.and_(hints.speciesid.in_(speciesid), hints.source == 'a2h', hints.seqnr == seqnr,
                                hints.start >= start, hints.end <= start + chunksize))
        exon_count = query.count()
        selected_intervals.append(maf_chunk)
        seen_exons += exon_count
        if seen_exons >= args.num_exons:
            break

    # run each chunk through augustus in training mode
    cgp_usage = tools.toilInterface.find_total_disk_usage([input_file_ids.fasta, input_file_ids.hints_db], buffer='4G')
    training_gffs = []
    for maf_chunk in selected_intervals:
        j = job.addChildJobFn(cgp, tree, maf_chunk, args, input_file_ids, training=True, memory='8G', disk=cgp_usage)
        training_gffs.append(j.rv())
    return job.addFollowOnJobFn(train_cgp, maf_chunks, tree, args, input_file_ids, training_gffs).rv()


def train_cgp(job, maf_chunks, tree, args, input_file_ids, training_gffs):
    """
    Trains CGP on a subset of alignments.
    :param maf_chunks: List of lists [chrom, start, chunksize, fileID]
    :param tree: fileID for tree from write_tree()
    :param args: input arguments
    :param input_file_ids: input file IDs
    :param training_gffs: List of fileIDs of chunks of CGP training
    :return: output of merge_results()
    """
    training_gff_files = [job.fileStore.readGlobalFile(x) for x in training_gffs]
    cmd = ['cat'] + training_gff_files
    combined = tools.fileOps.get_tmp_toil_file()
    tools.procOps.run_proc(cmd, stdout=combined)

    # run the training process
    params = tools.fileOps.get_tmp_toil_file()
    cmd = ['augustus',
           '--species={}'.format(args.species),
           '--treefile={}'.format(job.fileStore.readGlobalFile(tree)),
           '--refSpecies={}'.format(args.ref_genome),
           '--referenceFile={}'.format(job.fileStore.readGlobalFile(input_file_ids.gtf)),
           '--trainFeatureFile={}'.format(combined),
           '--param_outfile={}'.format(params)]
    tools.procOps.run_proc(cmd)
    input_file_ids.cgp_param = job.fileStore.writeGlobalFile(params)
    return job.addFollowOnJobFn(cgp_wrapper, maf_chunks, tree, args, input_file_ids).rv()


def cgp_wrapper(job, maf_chunks, tree, args, input_file_ids):
    """
    Wrapper for launching CGP jobs.
    :param maf_chunks: List of lists [chrom, start, chunksize, fileID]
    :param tree: fileID for tree from write_tree()
    :param args: input arguments
    :param input_file_ids: input file IDs
    :return: output of merge_results()
    """
    # results holds Promise objects
    # each Promise object will resolve to a tuple of gff_chunk_dict, stdout_file_id
    # cgp_job.rv():  key: genome, value: file handle to gff
    results = []
    cgp_usage = tools.toilInterface.find_total_disk_usage([input_file_ids.fasta, input_file_ids.hints_db], buffer='4G')

    for chrom, start, chunksize, maf_chunk in maf_chunks:
        # run AugustusCGP on alignment chunk
        cgp_job = job.addChildJobFn(cgp, tree, maf_chunk, args, input_file_ids, memory='8G', disk=cgp_usage)
        results.append([chrom, start, chunksize, cgp_job.rv()])

    # merge all gff files for alignment chunks to one gff for each species
    # results is a 3-member tuple of a joined genes list, a stdout file id dict and a file ID for the trained parameters
    # stdout_file_id dict is keyed by (chromosome, start, chunksize) tuples
    # for the joined genes dict its a dict keyed by genome and values are a 3 member tuple of:
    # [raw_gtf_file_id, joined_gtf_file_id, joined_gp_file_id]
    results = job.addFollowOnJobFn(merge_results, results, input_file_ids, memory='8G', disk='8G').rv()
    return results


def hal2maf(job, input_file_ids, genomes, ref_genome, annotate_ancestors, chrom, start, chunk_size):
    """
    exports hal to maf on a genomic region specified by (genome, seq, start, len)
    """
    hal = job.fileStore.readGlobalFile(input_file_ids.hal)
    maf_chunk = tools.fileOps.get_tmp_toil_file()
    genomes = ','.join(genomes)
    cmd = ['hal2maf', '--onlyOrthologs', '--refGenome', ref_genome, '--targetGenomes', genomes,
           '--refSequence', chrom, '--start', str(start), '--length', str(chunk_size), hal, maf_chunk]
    if not annotate_ancestors:
        cmd.append('--noAncestors')
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(maf_chunk)


def cgp(job, tree, maf_chunk, args, input_file_ids, training=False):
    """
    core function that runs AugustusCGP on one alignment chunk
    """
    genome_fofn = write_genome_fofn(job, input_file_ids.fasta)
    cgp_cfg = job.fileStore.readGlobalFile(input_file_ids.cgp_cfg)
    stdout = tools.fileOps.get_tmp_toil_file()

    cmd = ['augustus', '--dbhints=1', '--allow_hinted_splicesites=atac',
           '--extrinsicCfgFile={}'.format(cgp_cfg),
           '--species={}'.format(args.species),
           '--treefile={}'.format(job.fileStore.readGlobalFile(tree)),
           '--alnfile={}'.format(job.fileStore.readGlobalFile(maf_chunk)),
           '--dbaccess={}'.format(job.fileStore.readGlobalFile(input_file_ids.hints_db)),
           '--speciesfilenames={}'.format(genome_fofn),
           '--softmasking=1',
           '--exoncands={}'.format(1 if training else 0),
           '--alternatives-from-evidence=0',
           '--/CompPred/logreg=on',
           '--printOEs={}'.format(1 if training else 0),
           '--/CompPred/outdir={}'.format(os.getcwd())]
    if training is False:
        cmd.append('--optCfgFile={}'.format(job.fileStore.readGlobalFile(input_file_ids.cgp_param)))
    else:
        cmd.append('--printSampled=true')
    tools.procOps.run_proc(cmd, stdout=stdout)
    if training is True:
        cmd = ['cat', os.path.abspath('{}.sampled_GFs.gff'.format(args.ref_genome)),
               os.path.abspath('exonCands.{}.gff3'.format(args.ref_genome)),
               os.path.abspath('orthoExons.{}.gff3'.format(args.ref_genome))]
        combined_file = tools.fileOps.get_tmp_toil_file()
        tools.procOps.run_proc(cmd, stdout=combined_file)
        return job.fileStore.writeGlobalFile(combined_file)
    else:
        stdout_file_id = job.fileStore.writeGlobalFile(stdout)
        return {genome: job.fileStore.writeGlobalFile(genome + '.cgp.gff') for genome in args.genomes}, stdout_file_id


def merge_results(job, results, input_file_ids):
    """
    Results is a list of lists in the form [chrom, start, chunksize, (gff_chunk_dict, stdout_file_id)]
    gff_chunk is a dict of {genome: gff_file_id}
    Merges the results using joinGenes.
    """
    # reshape results into a dict of dicts:
    # {genome: (chrom, start, chunksize): gff_file_id
    gff_chunks_by_genome = collections.defaultdict(dict)
    stdout_file_ids = {}
    for chrom, start, chunksize, (gff_chunks, stdout_file_id) in results:
        stdout_file_ids[(chrom, start, chunksize)] = stdout_file_id
        for genome, gff_file_id in gff_chunks.iteritems():
            gff_chunks_by_genome[genome][(chrom, start, chunksize)] = gff_file_id
    results = {}
    for genome in gff_chunks_by_genome:
        j = job.addChildJobFn(join_genes, gff_chunks_by_genome[genome], memory='8G', disk='8G')
        results[genome] = j.rv()
    return results, stdout_file_ids, input_file_ids.cgp_param


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
    useful_lines = 0
    with open(raw_gtf_file, 'w') as raw_handle, open(raw_gtf_fofn, 'w') as fofn_handle:
        for (chrom, start, chunksize), chunk in gff_chunks.iteritems():
            local_path = job.fileStore.readGlobalFile(chunk)
            if os.environ.get('CAT_BINARY_MODE') == 'singularity':
                local_path = tools.procOps.singularify_arg(local_path)
            fofn_handle.write(local_path + '\n')
            raw_handle.write('## BEGIN CHUNK chrom: {} start: {} chunksize: {}\n'.format(chrom, start, chunksize))
            for line in open(local_path):
                if not line.startswith('#'):
                    useful_lines += 1
                raw_handle.write(line)

    # make sure CGP didn't fail entirely
    if useful_lines == 0:
        raise Exception('After running AugustusCGP, no gene predictions were made. Did you set `--augustus-species` '
                        'to a species with a trained model similar to your reference species? Please consult the '
                        'AUGUSTUS manual for more about the species flag.')

    join_genes_file = tools.fileOps.get_tmp_toil_file()
    join_genes_gp = tools.fileOps.get_tmp_toil_file()
    cmd = [['joingenes', '-f', raw_gtf_fofn, '-o', '/dev/stdout'],
           ['grep', '-P', '\tAUGUSTUS\t(exon|CDS|start_codon|stop_codon|tts|tss)\t'],
           ['sed', ' s/jg/augCGP-/g']]
    tools.procOps.run_proc(cmd, stdout=join_genes_file)

    # passing the joingenes output through gtfToGenePred then genePredToGtf fixes the sort order for homGeneMapping
    cmd = ['gtfToGenePred', '-genePredExt', join_genes_file, join_genes_gp]
    tools.procOps.run_proc(cmd)
    cmd = ['genePredToGtf', 'file', join_genes_gp, '-utr', '-honorCdsStat', '-source=augustusCGP', join_genes_file]
    tools.procOps.run_proc(cmd)

    joined_gtf_file_id = job.fileStore.writeGlobalFile(join_genes_file)
    raw_gtf_file_id = job.fileStore.writeGlobalFile(raw_gtf_file)
    joined_gp_file_id = job.fileStore.writeGlobalFile(join_genes_gp)
    return raw_gtf_file_id, joined_gtf_file_id, joined_gp_file_id


###
# Accessory functions
###


def write_tree(job,input_file_ids):
    """
    writes a file with the phylogenetic tree in NEWICK format
    """
    hal = job.fileStore.readGlobalFile(input_file_ids.hal)
    cmd = ['halStats', '--tree', hal]
    tree = tools.fileOps.get_tmp_toil_file()
    tools.procOps.run_proc(cmd, stdout=tree)
    return job.fileStore.writeGlobalFile(tree)


def write_genome_fofn(job, fasta_file_ids):
    """
    writes a file with the location of the fasta files, e.g.

    galGal4 /path/to/genome/galGal4.fa
    hg38    /path/to/genome/hg38.fa
    mm10    /path/to/genome/mm10.fa
    rn6     /path/to/genome/rn6.fa
    ...

    These files are loaded from the fileStore
    """
    genome_fofn = tools.fileOps.get_tmp_toil_file()
    with open(genome_fofn, 'w') as outf:
        for genome, file_id in fasta_file_ids.iteritems():
            local_path = job.fileStore.readGlobalFile(file_id)
            if os.environ.get('CAT_BINARY_MODE') == 'singularity':
                local_path = tools.procOps.singularify_arg(local_path)
            tools.fileOps.print_row(outf, [genome, local_path])
    return genome_fofn


