"""
Toil pipeline to align all transcripts in a source genePred to a target genePred. For non-CGP transcripts,
this can be determined by using tools.nameConversions.strip_alignment_numbers() on the name field. For CGP transcripts,
which have new IDs, we use the name2 field which will have assigned a gene ID to try and align to all protein coding
transcripts associated with that gene ID.

Alignment is performed in a few different ways:
1. For each CGP transcript, the in-frame CDS will be aligned using PRANK to the in-frame CDS of each protein-coding
transcript of the assigned parental gene.
2. For each transMap transcript, it will be aligned via MUSCLE to the assigned parent. If the parent is protein coding
then the transcript will also undergo in-frame CDS alignment via PRANK.
3. For each AugustusTM(R) transcript, it will be aligned both via MUSCLE and PRANK.

"""
import logging
import collections
import os
import itertools

from toil.job import Job
from toil.common import Toil

import tools.bio
import tools.psl
import tools.dataOps
import tools.fileOps
import tools.procOps
import tools.transcripts
import tools.toilInterface
import tools.sqlInterface
import tools.nameConversions


def align_transcripts(args, toil_options):
    """
    Main entry function for transcript alignment toil pipeline
    :param args: dictionary of arguments from CAT
    :param toil_options: toil options Namespace object
    """
    with Toil(toil_options) as toil:
        if not toil.options.restart:
            # assume that both fasta's have been flattened
            tx_file_ids = tools.toilInterface.write_fasta_to_filestore(toil, args['ref_genome_fasta'])
            ref_genome_fasta_file_id, ref_genome_fasta_gdx_file_id, ref_genome_fasta_flat_file_id = tx_file_ids
            genome_file_ids = tools.toilInterface.write_fasta_to_filestore(toil, args['genome_fasta'])
            genome_fasta_file_id, genome_fasta_gdx_file_id, genome_fasta_flat_file_id = genome_file_ids
            input_file_ids = {'ref_genome_fasta': ref_genome_fasta_file_id,
                              'ref_genome_gdx': ref_genome_fasta_gdx_file_id,
                              'ref_genome_flat': ref_genome_fasta_flat_file_id,
                              'genome_fasta': genome_fasta_file_id,
                              'genome_gdx': genome_fasta_gdx_file_id,
                              'genome_flat': genome_fasta_flat_file_id,
                              'annotation_gp': toil.importFile('file://' + args['annotation_gp']),
                              'annotation_db': toil.importFile('file://' + args['annotation_db']),
                              'modes': {}}
            for mode in args['modes']:
                input_file_ids['modes'][mode] = toil.importFile('file://' + args['modes'][mode]['gp'])
            job = Job.wrapJobFn(setup, args, input_file_ids)
            results_file_ids = toil.start(job)
        else:
            results_file_ids = toil.restart()
        for file_path, file_id in results_file_ids.iteritems():
            tools.fileOps.ensure_file_dir(file_path)
            toil.exportFile(file_id, 'file://' + file_path)


def setup(job, args, input_file_ids):
    """
    First function for align_transcripts pipeline. Splits up the genePred entries into chunks that will be aligned
    with BLAT.
    :param args: dictionary of arguments from CAT
    :param input_file_ids: dictionary of fileStore file IDs for the inputs to this pipeline
    """
    job.fileStore.logToMaster('Beginning Align Transcripts run on {}'.format(args['genome']), level=logging.INFO)
    chunk_size = 100
    cgp_chunk_size = 20  # cgp has multiple alignments
    # load all fileStore files necessary
    annotation_gp = job.fileStore.readGlobalFile(input_file_ids['annotation_gp'])
    annotation_db = job.fileStore.readGlobalFile(input_file_ids['annotation_db'])
    # we have to explicitly place fasta, flat file and gdx with the correct naming scheme for pyfasta
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids['genome_fasta'],
                                                                 input_file_ids['genome_gdx'],
                                                                 input_file_ids['genome_flat'],
                                                                 prefix='genome', upper=False)
    ref_genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids['ref_genome_fasta'],
                                                                     input_file_ids['ref_genome_gdx'],
                                                                     input_file_ids['ref_genome_flat'],
                                                                     prefix='ref_genome', upper=False)
    # load required reference data into memory
    tx_biotype_map = tools.sqlInterface.get_transcript_biotype_map(annotation_db, args['ref_genome'])
    ref_transcript_dict = tools.transcripts.get_gene_pred_dict(annotation_gp)
    # will hold a mapping of output file paths to lists of Promise objects containing output
    results = collections.defaultdict(list)
    # start generating chunks of the transMap/Augustus genePreds which we know the 1-1 alignment for
    for mode in ['transMap', 'augTM', 'augTMR']:
        if mode not in args['modes']:
            continue
        # output file paths
        muscle_path = args['modes'][mode]['MUSCLE']
        prank_path = args['modes'][mode]['PRANK']
        # begin loading transcripts and sequences
        gp_path = job.fileStore.readGlobalFile(input_file_ids['modes'][mode])
        transcript_dict = tools.transcripts.get_gene_pred_dict(gp_path)
        muscle_seq_iter = get_muscle_sequences(transcript_dict, ref_transcript_dict, genome_fasta, ref_genome_fasta)
        for chunk in tools.dataOps.grouper(muscle_seq_iter, chunk_size):
            j = job.addChildJobFn(run_muscle_chunk, chunk)
            results[muscle_path].append(j.rv())
        prank_seq_iter = get_prank_sequences(transcript_dict, ref_transcript_dict, genome_fasta, ref_genome_fasta,
                                             tx_biotype_map)
        for chunk in tools.dataOps.grouper(prank_seq_iter, chunk_size):
            j = job.addChildJobFn(run_prank_chunk, chunk)
            results[prank_path].append(j.rv())
    if 'augCGP' in args['modes']:
        cgp_prank_path = args['modes']['augCGP']['PRANK']
        # CGP transcripts have multiple assignments based on the name2 identifier, which contains a gene ID
        gene_tx_map = tools.sqlInterface.get_gene_transcript_map(annotation_db, args['ref_genome'])
        tx_biotype_map = tools.sqlInterface.get_transcript_biotype_map(annotation_db, args['ref_genome'])
        augustus_cgp_gp = job.fileStore.readGlobalFile(input_file_ids['augustus_cgp_gp'])
        cgp_transcript_dict = tools.transcripts.get_gene_pred_dict(augustus_cgp_gp)
        cgp_transcript_seq_iter = get_cgp_sequences(cgp_transcript_dict, ref_transcript_dict, genome_fasta,
                                                    ref_genome_fasta, gene_tx_map, tx_biotype_map)
        for chunk in tools.dataOps.grouper(cgp_transcript_seq_iter, cgp_chunk_size):
            j = job.addChildJobFn(run_prank_chunk, chunk)
            results[cgp_prank_path].append(j.rv())
    if len(results) == 0:
        err_msg = 'Align Transcripts pipeline did not detect any input genePreds for {}'.format(args['genome'])
        raise RuntimeError(err_msg)
    # convert the results Promises into resolved values
    return job.addFollowOnJobFn(merge, results, args).rv()


def get_muscle_sequences(transcript_dict, ref_transcript_dict, genome_fasta, ref_genome_fasta):
    """Generator that yields a tuple of (tx_id, tx_seq, ref_tx_id, ref_tx_seq)"""
    for tx_id, tx in transcript_dict.iteritems():
        ref_tx_id = tools.nameConversions.strip_alignment_numbers(tx_id)
        ref_tx = ref_transcript_dict[ref_tx_id]
        tx_seq = tx.get_mrna(genome_fasta)
        ref_tx_seq = ref_tx.get_mrna(ref_genome_fasta)
        yield tx_id, tx_seq, ref_tx_id, ref_tx_seq


def get_prank_sequences(transcript_dict, ref_transcript_dict, genome_fasta, ref_genome_fasta, tx_biotype_map):
    """Generator that yields a tuple of (tx_id, tx_seq, ref_tx_id, ref_tx_seq)"""
    for tx_id, tx in transcript_dict.iteritems():
        ref_tx_id = tools.nameConversions.strip_alignment_numbers(tx_id)
        ref_tx = ref_transcript_dict[ref_tx_id]
        biotype = tx_biotype_map[ref_tx_id]
        if biotype == 'protein_coding':
            tx_seq = tx.get_cds(genome_fasta, in_frame=True)
            ref_tx_seq = ref_tx.get_cds(ref_genome_fasta, in_frame=True)
            assert len(tx_seq) % 3 == 0, tx_id
            assert len(ref_tx_seq) % 3 == 0, ref_tx_id
            # don't let PRANK try to align very short sequences
            if len(ref_tx_seq) > 50 and len(tx_seq) > 50:
                yield tx_id, tx_seq, ref_tx_id, ref_tx_seq


def get_cgp_sequences(transcript_dict, ref_transcript_dict, genome_fasta, ref_genome_fasta, gene_tx_map,
                      tx_biotype_map):
    """
    Generator for CGP transcripts. Same as get_prank_sequences, but will resolve name2 field into all target transcripts
    """
    for cgp_id, tx in transcript_dict.iteritems():
        if 'jg' in tx.name2:
            continue  # this transcript was not assigned any parents
        ref_tx_ids = gene_tx_map[tx.name2]
        tx_seq = tx.get_cds(genome_fasta, in_frame=True)
        assert len(tx_seq) % 3 == 0, cgp_id
        for ref_tx_id in ref_tx_ids:
            biotype = tx_biotype_map[ref_tx_id]
            if biotype != 'protein_coding':
                continue
            ref_tx = ref_transcript_dict[ref_tx_id]
            ref_tx_seq = ref_tx.get_cds(ref_genome_fasta, in_frame=True)
            assert len(ref_tx_seq) % 3 == 0, ref_tx_id
            if len(ref_tx_seq) > 50 and len(tx_seq) > 50:
                yield cgp_id, tx_seq, ref_tx_id, ref_tx_seq


def run_prank_chunk(job, chunk):
    """
    Runs an alignment chunk through PRANK for coding cds alignment. Requires post-processing because PRANK handles
    stop codons in dumb ways.
    :param chunk: List of (tx_id, tx_seq, ref_tx_id, ref_tx_seq) tuples
    :return: List of PSL output
    """
    tmp_fasta = tools.fileOps.get_tmp_toil_file()
    results = []
    cmd = ['prank', '-codon', '-DNA', '-d={}'.format(tmp_fasta)]
    prank_out_file = 'output.best.fas'
    for tx_id, tx_seq, ref_tx_id, ref_tx_seq in chunk:
        with open(tmp_fasta, 'w') as outf:
            tools.bio.write_fasta(outf, ref_tx_id, ref_tx_seq)
            tools.bio.write_fasta(outf, tx_id, tx_seq)
        tools.procOps.run_proc(cmd)
        # prank may fail to find any alignments
        if not os.path.exists(prank_out_file):
            continue
        r = parse_prank_output(prank_out_file, tx_seq, ref_tx_seq)
        results.append(r)
        os.remove(prank_out_file)  # we need to remove this file in case the next iteration fails
    return results


def run_muscle_chunk(job, chunk):
    """
    Runs an alignment chunk through muscle for non-coding or coding mRNA sequences
    :param chunk: List of (tx_id, tx_seq, ref_tx_id, ref_tx_seq) tuples
    :return: List of PSL output
    """
    tmp_fasta = tools.fileOps.get_tmp_toil_file()
    results = []
    cmd = ['muscle', '-in', tmp_fasta]
    for tx_id, tx_seq, ref_tx_id, ref_tx_seq in chunk:
        with open(tmp_fasta, 'w') as outf:
            tools.bio.write_fasta(outf, ref_tx_id, ref_tx_seq)
            tools.bio.write_fasta(outf, tx_id, tx_seq)
        r = tools.procOps.call_proc(cmd, keepLastNewLine=True)
        results.append(r)
    return results


def merge(job, results, args):
    """
    Merge together chain files.
    :param results: dict of list of promises from each alignment chunk for each category
    :param args: arguments to the pipeline
    :return:
    """
    job.fileStore.logToMaster('Merging Alignment output for {}'.format(args['genome']), level=logging.INFO)
    results_file_ids = {}
    for gp_category, result_list in results.iteritems():
        tmp_results_file = tools.fileOps.get_tmp_toil_file(suffix='gz')
        with tools.fileOps.opengz(tmp_results_file, 'w') as outf:
            for line in itertools.chain.from_iterable(result_list):  # results is list of lists
                outf.write(line + '\n')
        results_file_ids[gp_category] = job.fileStore.writeGlobalFile(tmp_results_file)
    return results_file_ids


###
# Helper functions
###


def parse_prank_output(prank_out_file, tx_seq, ref_tx_seq):
    """
    Converts the output from PRANK into a sane alignment. PRANK likes to mess up stop codons, like so:

    Example 1. Pruning an in-frame terminal stop codon from the reference:

    Input:
    >ref
    ATGGGGTGA
    >tgt
    ATGGGGTGC

    Output:
    >ref
    ATGGGG---
    >tgt
    ATGGGGTGC

    Desired output:
    >ref
    ATGGGGTGA
    >tgt
    ATGGGGTGC


    Example 2: Does not prune an in-frame stop codon if it exists on both:
    Input, Output, Desired output:
    >ref
    ATGGGGTGA
    >tgt
    ATGGGGTGA


    Example 3: Adds Ns when an in-frame stop codon is aligned over:
    Input, desired output:
    >ref
    ATGGGGTGATAC
    >tgt
    ATGGGGTGATAC

    Output:
    >ref
    ATGGGGNNNTAC
    >tgt
    ATGGGGNNNTAC


    Example 4: The same as example 3, but where the stop only exists on one (more common):
    Input, desired output:
    >ref
    ATGGGGTGATAC
    >tgt
    ATGGGGTCATAC

    Output:
    >ref
    ATGGGGNNNTAC
    >tgt
    ATGGGGTCATAC

    :param prank_out_file: Path to PRANK output. Generally 'output.best.fas'
    :param tx_seq: Sequence of the target transcript for comparison to PRANK output
    :param ref_tx_seq: Sequence of the reference transcript for comparison to PRANK output
    :return: Munged string mimicking PRANK output after analysis
    """
    ref_name, ref_aln, tgt_name, tgt_aln = list(itertools.chain.from_iterable(tools.bio.read_fasta(prank_out_file)))
    # fasta naming
    ref_name = '>' + ref_name
    tgt_name = '>' + tgt_name
    # first, we test to see if PRANK broke anything
    if validate_prank(ref_aln, ref_tx_seq) and validate_prank(tgt_aln, tx_seq):
        r = '\n'.join([ref_name, ''.join(ref_aln), tgt_name, ''.join(tgt_aln), ''])
        return r
    # start walking the alignment, inserting the correct base where it should be
    new_ref_aln = []
    new_tgt_aln = []
    ref_pos = 0
    tgt_pos = 0
    for aln_pos, (ref_aln_base, tgt_aln_base) in enumerate(itertools.izip(*[ref_aln, tgt_aln])):
        if ref_aln_base != '-':
            if ref_aln_base != ref_tx_seq[ref_pos]:
                new_ref_aln.append(ref_tx_seq[ref_pos])
            else:
                new_ref_aln.append(ref_aln_base)
            ref_pos += 1
        else:
            new_ref_aln.append('-')
        if tgt_aln_base != '-':
            if tgt_aln_base != tx_seq[tgt_pos]:
                new_tgt_aln.append(tx_seq[tgt_pos])
            else:
                new_tgt_aln.append(tgt_aln_base)
            tgt_pos += 1
        else:
            new_tgt_aln.append('-')
    if ref_pos != len(ref_tx_seq):
        # we need to add back in the stop codon
        assert len(ref_tx_seq) - ref_pos == 3
        new_ref_aln.extend(ref_tx_seq[ref_pos:])
    if tgt_pos != len(tx_seq):
        assert len(tx_seq) - tgt_pos == 3
        new_tgt_aln.extend(tx_seq[tgt_pos:])
    # we need to add in gaps for the shorter one
    if len(new_ref_aln) < len(new_tgt_aln):
        new_ref_aln.extend(['-'] * (len(new_tgt_aln) - len(new_ref_aln)))
    if len(new_tgt_aln) < len(new_ref_aln):
        new_tgt_aln.extend(['-'] * (len(new_ref_aln) - len(new_tgt_aln)))
    new_ref_aln = ''.join(new_ref_aln)
    new_tgt_aln = ''.join(new_tgt_aln)
    assert new_ref_aln.replace('-', '') == ref_tx_seq, (new_ref_aln, ref_tx_seq, new_tgt_aln, tx_seq)
    assert new_tgt_aln.replace('-', '') == tx_seq, (new_ref_aln, ref_tx_seq, new_tgt_aln, tx_seq)
    assert len(new_ref_aln) == len(new_tgt_aln)
    r = '\n'.join([ref_name, new_ref_aln, tgt_name, new_tgt_aln, ''])
    return r


def validate_prank(aln, seq, rm='-N'):
    """
    Uses string.translate to remove the dashes and Ns from a pair of sequences and asks if they are the same
    :param aln: String derived from the output of PRANK
    :param seq: String derived from the input to PRANK
    :param rm: characters to remove. Defaults to '-N'
    :return: boolean (true if things are messed up)
    """
    return aln.translate(None, rm) == seq.translate(None, rm)
