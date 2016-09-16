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

from toil.job import Job
from toil.common import Toil

import logging
import argparse
import itertools
import os
import collections

import tools.fileOps
import tools.dataOps
import tools.procOps
import tools.intervals
import tools.transcripts
import tools.sqlInterface


###
# AugustusCGP pipeline section
###


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
            input_file_ids.cgp_param = toil.importFile('file://' + args.augustus_cgp_param)
            input_file_ids.ref_genome_db = toil.importFile('file://' + args.ref_genome_db)
            input_file_ids.fasta = {genome: toil.importFile('file://' + fasta)
                                    for genome, fasta in args.fasta_files.iteritems()}
            input_file_ids.tm_gps = {genome: toil.importFile('file://' + tm_gp)
                                     for genome, tm_gp in args.tm_gps.iteritems()}
            if args.augustus_cgp_cfg is not None:
                cgp_cfg_file_id = toil.importFile('file://' + args.augustus_cgp_cfg)
                input_file_ids.cgp_cfg = cgp_cfg_file_id
            job = Job.wrapJobFn(setup, args, input_file_ids, memory='8G')
            results = toil.start(job)
        else:
            results = toil.restart()
        for genome in results:
            tools.fileOps.ensure_file_dir(args.augustus_cgp_gtf[genome])
            toil.exportFile(results[genome], 'file://' + args.augustus_cgp_gtf[genome])


def setup(job, args, input_file_ids):
    """
    Entry function for running AugustusCGP.
    HAL alignment is converted to MAF format and splitted into overlapping
    alignment chunks for parallel computing. Each alignment chunk is one child process.
    Gene predictions on alignment chunks are subsequently merged into one gff for each species.
    For merging of the gene sets, the auxiliary tool 'joingenes' from the Augustus tool package is used.
    """
    job.fileStore.logToMaster('Preparing input files for AugustusCGP', level=logging.INFO)

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
    for i, (chrom, chromSize) in enumerate(tools.fileOps.iter_lines(chromSizes)):
        start = 0
        while start + chunkSize < int(chromSize):
            aliChunks.append([chrom, start, chunkSize])
            start = start + chunkSize - overlap                  # start of next alignment chunk
        aliChunks.append([chrom, start, int(chromSize) - start])   # last alignment chunk

    for chrom, start, end in aliChunks:
        # string "genome.chrom:start-end"
        genomic_region = '{}.{}:{}-{}'.format(args.ref_genome, chrom, start, start + end - 1)
        # export alignment chunks from hal to maf
        j = job.addChildJobFn(hal2maf, input_file_ids, args.ref_genome, chrom, start, end, genomic_region,
                              memory='8G')
        mafChunk = j.rv()
        # run AugustusCGP on alignment chunk
        cgp_job = j.addFollowOnJobFn(cgp, tree, mafChunk, args, input_file_ids, genomic_region, memory='8G')
        gffChunk = cgp_job.rv()
        gffChunks.append(gffChunk)

    # merge all gff files for alignment chunks to one gff for each species
    mergedGffs = job.addFollowOnJobFn(merge_results, args, input_file_ids, gffChunks, memory='8G').rv()
    return mergedGffs


def hal2maf(job, input_file_ids, refGenome, chrom, start, chunkSize, genomic_region):
    """
    exports hal to maf on a genomic region specified by (genome, seq, start, len)
    """
    job.fileStore.logToMaster('Running hal2maf on {}'.format(genomic_region), level=logging.INFO)
    hal = job.fileStore.readGlobalFile(input_file_ids.hal)
    mafChunk = tools.fileOps.get_tmp_toil_file()
    cmd = ['hal2maf', '--noAncestors', '--noDupes', '--refGenome', refGenome,
           '--refSequence', chrom, '--start', start, '--length', chunkSize, hal, mafChunk]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(mafChunk)


def cgp(job, tree, mafChunk, args, input_file_ids, genomic_region):
    """
    core function that runs AugustusCGP on one alignment chunk
    """
    job.fileStore.logToMaster('Running AugustusCGP on {}'.format(genomic_region), level=logging.INFO)

    genomeFofn = writeGenomeFofn(job, input_file_ids.fasta)

    opt_param = []  # optional AugustusCGP parameters
    if args.augustus_cgp_cfg is not None:
        # extrinsic config file
        # for now, let's assume RNA-Seq evidence whenever an extrinsic
        # config file is given, e.g. by default we turn on
        # --UTR=1
        # --allow_hinted_splicesites=atac
        cgp_cfg = job.fileStore.readGlobalFile(args.cgp_cfg)
        opt_param += ['--extrinsicCfgFile={}'.format(cgp_cfg),
                      '--dbhints=1', '--UTR=1', '--allow_hinted_splicesites=atac']

    cmd = ['augustus'] + opt_param + ['--species={}'.format(args.species),
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
    mergedGffs = {}
    for genome in args.genomes:
        # merge all gffChunks of one genome
        genome_gffChunks = [d[genome] for d in gffChunks]
        j = job.addChildJobFn(joinGenes, genome, args, input_file_ids, genome_gffChunks, memory='8G')
        mergedGffs[genome] = j.rv()
    return mergedGffs


def joinGenes(job, genome, args, input_file_ids, gffChunks):
    """
    uses the auxiliary tool 'joingenes' from the
    Augustus package to intelligently merge gene sets
    - removes duplicated Txs or truncated Txs that are contained in other Txs (trivial)
    - fixes truncated Txs at alignment boundaries,
      e.g. by merging them with other Txs (non trivial, introduces new Txs)
    
    Calls out to the parental gene assignment pipeline
    """
    job.fileStore.logToMaster('Merging GFFs for {}'.format(genome), level=logging.INFO)

    fofn = tools.fileOps.get_tmp_toil_file()
    with open(fofn, 'w') as outf:
        for chunk in gffChunks:
            local_path = job.fileStore.readGlobalFile(chunk)
            outf.write(local_path + '\n')

    jg = tools.fileOps.get_tmp_toil_file()
    cmd = [['joingenes', '-f', fofn, '-o', '/dev/stdout'],
       ['grep', '-P', '\tAUGUSTUS\t(exon|CDS|start_codon|stop_codon|tts|tss)\t']]
    tools.procOps.run_proc(cmd, stdout=jg)
    joined_file_id = job.fileStore.writeGlobalFile(jg)
    j = job.addFollowOnJobFn(assign_parents, args, genome, input_file_ids, joined_file_id, memory='8G')
    return j.rv()


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


###
# Parent gene assignment section
###


def assign_parents(job, args, genome, input_file_ids, joined_gff_file_id):
    """
    Main function for assigning parental genes. Parental gene assignment methodology:
    A) Each CGP transcript is evaluated for overlapping any transMap transcripts. Overlap is defined as having at least
    1 exonic base shared.
    B) If a CGP transcript is assigned to more than one gene, then this is attempted to be resolved by first:
        1) Determine if only one of the possible genes is coding. If so, make the assignment.
        2) Use the Jaccard distance metric to attempt to find one best hit.
        If none of the above are satisfied, the CGP transcript is discarded as unresolvable.

    TODO: this process operates on genePreds, and so I needlessly transform to genePred then back. This could introduce
    errors. The better way would be a proper GFF parser. Mark's Gff3 parser does not work.
    """
    job.fileStore.logToMaster('Assigning parental genes for {}'.format(genome), level=logging.INFO)
    ref_genome_db = job.fileStore.readGlobalFile(input_file_ids.ref_genome_db)
    gene_biotype_map = tools.sqlInterface.get_gene_biotype_map(ref_genome_db)
    tm_gp_file = job.fileStore.readGlobalFile(input_file_ids.tm_gps[genome])
    transmap_dict = tools.transcripts.get_gene_pred_dict(tm_gp_file)

    # convert GFF to genePred
    cgp_gp = tools.fileOps.get_tmp_toil_file()
    joined_gff = job.fileStore.readGlobalFile(joined_gff_file_id)
    cmd = ['gtfToGenePred', '-genePredExt', joined_gff, cgp_gp]
    tools.procOps.run_proc(cmd)
    cgp_dict = tools.transcripts.get_gene_pred_dict(cgp_gp)
    tm_chrom_dict = create_chrom_dict(transmap_dict)
    cgp_chrom_dict = create_chrom_dict(cgp_dict)

    final_gps = []
    for chrom, tm_tx_by_chromosome in tm_chrom_dict.iteritems():
        for cgp_chunk in tools.dataOps.grouper(cgp_chrom_dict[chrom].iteritems(), 70):
            j = job.addChildJobFn(assign_parent_chunk, tm_tx_by_chromosome, cgp_chunk, gene_biotype_map)
            final_gps.append(j.rv())
    return job.addFollowOnJobFn(merge_parent_assignment_chunks, final_gps).rv()


def assign_parent_chunk(job, tm_tx_by_chromosome, cgp_chunk, gene_biotype_map):
    """
    Runs a chunk of CGP transcripts on the same chromosome as all transMap transcripts in tm_tx_by_chromosome
    :param tm_tx_by_chromosome: dictionary of GenePredTranscript objects all on the same chromosome
    :param cgp_chunk: Iterable of (cgp_tx_id, cgp_tx) tuples to be analyzed
    :param gene_biotype_map: dictionary mapping gene IDs to biotype
    :return: list of GenePredTranscript objects which have been resolved
    """
    resolved_txs = []
    for cgp_tx_id, cgp_tx in cgp_chunk:
        overlapping_tm_txs = find_tm_overlaps(cgp_tx, tm_tx_by_chromosome)
        gene_ids = {tx.name2 for tx in overlapping_tm_txs}
        if len(gene_ids) == 0:
            gene_name = cgp_tx.name.split('.')[0]
        elif len(gene_ids) == 1:
            gene_name = list(gene_ids)[0]
        else:
            gene_name = resolve_multiple_genes(cgp_tx, overlapping_tm_txs, gene_biotype_map)
        if gene_name is not None:  # we can resolve this transcript
            cgp_tx.name2 = gene_name
            resolved_txs.append(cgp_tx)
    return resolved_txs


def merge_parent_assignment_chunks(job, final_gps):
    """
    Merge the chunks of transcripts produced by assign_parent_chunk, converting back to GFF
    :param final_gps: list of lists of GenePredTranscript objects
    :return: fileStore ID to a output GFF
    """
    out_gp = tools.fileOps.get_tmp_toil_file()
    out_gff = tools.fileOps.get_tmp_toil_file()
    with open(out_gp, 'w') as outf:
        for rec in itertools.chain.from_iterable(final_gps):
            tools.fileOps.print_row(outf, rec.get_gene_pred())
    cmd = ['genePredToGtf', '-utr', '-honorCdsStat', '-source=AugustusCGP', 'file', out_gp, out_gff]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(out_gff)


###
# Parent gene assignment functions
###


def create_chrom_dict(tx_dict):
    """
    Split up a dictionary of Transcript objects by chromosome
    """
    chrom_dict = collections.defaultdict(dict)
    for tx_id, tx in tx_dict.iteritems():
        chrom_dict[tx.chromosome][tx_id] = tx
    return chrom_dict


def find_tm_overlaps(cgp_tx, tm_tx_dict):
    """Find overlap with transMap transcripts first on a genomic scale then an exonic scale"""
    r = []
    for tx in tm_tx_dict.itervalues():
        if tx.interval.intersection(cgp_tx.interval) is not None:
            # make sure that we have exon overlap
            if ensure_exon_overlap(tx, cgp_tx) is True:
                r.append(tx)
    return r


def ensure_exon_overlap(tx, cgp_tx):
    """Do these two transcripts have at least 1 exonic base of overlap?"""
    for tm_exon in tx.exon_intervals:
        for cgp_exon in cgp_tx.exon_intervals:
            if tm_exon.overlap(cgp_exon) is True:
                return True
    return False


def resolve_multiple_genes(cgp_tx, overlapping_tm_txs, gene_biotype_map, min_distance=0.2):
    """
    Resolve multiple assignments based on the following rules:
    If a CGP has multiple assignments, see if only 1 is protein_coding. If so, assign it.
    If the difference in Jaccard scores is >=min_distance, assign it to the higher score.
    If neither of these are satisfiable, discard the transcript.
    """
    filtered = [tm_tx for tm_tx in overlapping_tm_txs if gene_biotype_map[tm_tx.name2] == 'protein_coding']
    if len(filtered) == 0:
        return None
    elif len(filtered) == 1:
        return filtered[0].name2
    # attempt the Jaccard metric
    scores = calculate_jaccard(cgp_tx, filtered)
    high_score = max(scores.itervalues())
    if all(high_score - x >= min_distance for x in scores.itervalues() if x != high_score):
        best = sorted(scores.iteritems(), key=lambda (gene_id, score): score)[-1]
        results = best[0]
    else:
        results = None
    return results


def calculate_jaccard(cgp_tx, filtered_overlapping_tm_txs):
    """calculates the Jaccard distance metric."""
    results = collections.defaultdict(float)
    with tools.fileOps.TemporaryFilePath() as cgp:
        with open(cgp, 'w') as outf:
            tools.fileOps.print_row(outf, cgp_tx.get_bed())
        with tools.fileOps.TemporaryFilePath() as tm:
            for tm_tx in filtered_overlapping_tm_txs:
                with open(tm, 'w') as outf:
                    tools.fileOps.print_row(outf, tm_tx.get_bed())
                cmd = ['bedtools', 'jaccard', '-s', '-a', cgp, '-b', tm]
                r = tools.procOps.call_proc_lines(cmd)
                j = float(r[-1].split()[-2])
                results[tm_tx.name2] = max(results[tm_tx.name2], j)
    return results
