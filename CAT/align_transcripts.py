"""
Toil pipeline to align all transcripts in a source genePred to a target genePred. For non-CGP transcripts,
this can be determined by using tools.nameConversions.strip_alignment_numbers() on the name field. For CGP transcripts,
which have new IDs, we use the name2 field which will have assigned a gene ID to try and align to all protein coding
transcripts associated with that gene ID.
"""
import logging
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
            tm_gp_file_id = toil.importFile('file://' + args['tm_gp'])
            annotation_gp_file_id = toil.importFile('file://' + args['annotation_gp'])
            annotation_db_file_id = toil.importFile('file://' + args['annotation_db'])
            input_file_ids = {'ref_genome_fasta': ref_genome_fasta_file_id,
                              'ref_genome_gdx': ref_genome_fasta_gdx_file_id,
                              'ref_genome_flat': ref_genome_fasta_flat_file_id,
                              'genome_fasta': genome_fasta_file_id,
                              'genome_gdx': genome_fasta_gdx_file_id,
                              'genome_flat': genome_fasta_flat_file_id,
                              'tm_gp': tm_gp_file_id,
                              'annotation_gp': annotation_gp_file_id,
                              'annotation_db': annotation_db_file_id}
            if 'augustus_gp' in args:
                augustus_gp_file_id = toil.importFile('file://' + args['augustus_gp'])
                input_file_ids['augustus_gp'] = augustus_gp_file_id
            if 'augustus_cgp_gp' in args:
                augustus_cgp_gp_file_id = toil.importFile('file://' + args['augustus_cgp_gp'])
                input_file_ids['augustus_cgp_gp'] = augustus_cgp_gp_file_id
            job = Job.wrapJobFn(setup, args, input_file_ids)
            results_file_id = toil.start(job)
        else:
            results_file_id = toil.restart()
        tools.fileOps.ensure_file_dir(args['alignment_psl'])
        toil.exportFile(results_file_id, 'file://' + args['alignment_psl'])


def setup(job, args, input_file_ids):
    """
    First function for align_transcripts pipeline. Splits up the genePred entries into chunks that will be aligned
    with BLAT.
    :param args: dictionary of arguments from CAT
    :param input_file_ids: dictionary of fileStore file IDs for the inputs to this pipeline
    """
    job.fileStore.logToMaster('Beginning AlignTranscripts run on {}'.format(args['genome']), level=logging.INFO)
    chunk_size = 100
    cgp_chunk_size = 20  # CGP will have multiple alignments per transcript
    # load all fileStore files necessary
    tm_gp = job.fileStore.readGlobalFile(input_file_ids['tm_gp'])
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
    # start generating chunks of the transMap/Augustus genePreds which we know the 1-1 alignment for
    tx_biotype_map = tools.sqlInterface.get_transcript_biotype_map(annotation_db, args['ref_genome'])
    results = []
    gp_file_handles = [open(tm_gp)]
    if 'augustus_gp' in args:
        augustus_gp = job.fileStore.readGlobalFile(input_file_ids['augustus_gp'])
        gp_file_handles.append(open(augustus_gp))
    gp_iter = itertools.chain.from_iterable(gp_file_handles)
    transcript_dict = tools.transcripts.get_gene_pred_dict(gp_iter)
    ref_transcript_dict = tools.transcripts.get_gene_pred_dict(annotation_gp)
    transcript_seq_iter = get_sequences(transcript_dict, ref_transcript_dict, genome_fasta, ref_genome_fasta,
                                        tx_biotype_map)
    for i, chunk in enumerate(tools.dataOps.grouper(transcript_seq_iter, chunk_size)):
        j = job.addChildJobFn(run_alignment_chunk, i, args, chunk)
        results.append(j.rv())
    if 'augustus_cgp_gp' in args:
        # CGP transcripts have multiple assignments based on the name2 identifier, which contains a gene ID
        gene_tx_map = tools.sqlInterface.get_gene_transcript_map(annotation_db, args['ref_genome'])
        tx_biotype_map = tools.sqlInterface.get_transcript_biotype_map(annotation_db, args['ref_genome'])
        augustus_cgp_gp = job.fileStore.readGlobalFile(input_file_ids['augustus_cgp_gp'])
        cgp_transcript_dict = tools.transcripts.get_gene_pred_dict(augustus_cgp_gp)
        cgp_transcript_seq_iter = get_cgp_sequences(cgp_transcript_dict, ref_transcript_dict, genome_fasta,
                                                    ref_genome_fasta, gene_tx_map, tx_biotype_map)
        for i, chunk in enumerate(tools.dataOps.grouper(cgp_transcript_seq_iter, cgp_chunk_size), start=i):
            j = job.addChildJobFn(run_alignment_chunk, i, args, chunk)
            results.append(j.rv())
    return job.addFollowOnJobFn(merge, results, args).rv()


def get_sequences(transcript_dict, ref_transcript_dict, genome_fasta, ref_genome_fasta, tx_biotype_map):
    """Generator that yields a tuple of (tx_id, tx_seq, ref_tx_id, ref_tx_seq)"""
    for tx_id, tx in transcript_dict.iteritems():
        ref_tx_id = tools.nameConversions.strip_alignment_numbers(tx_id)
        ref_tx = ref_transcript_dict[ref_tx_id]
        if tx_biotype_map[ref_tx_id] == 'protein_coding':
            tx_seq = tx.get_cds(genome_fasta)
            ref_tx_seq = ref_tx.get_cds(ref_genome_fasta)
        else:
            tx_seq = tx.get_mrna(genome_fasta)
            ref_tx_seq = ref_tx.get_mrna(ref_genome_fasta)
        yield tx_id, tx_seq, ref_tx_id, ref_tx_seq


def get_cgp_sequences(transcript_dict, ref_transcript_dict, genome_fasta, ref_genome_fasta, gene_tx_map,
                      tx_biotype_map):
    """Generator for CGP transcripts. Same as get_sequences, but will resolve name2 field into all target transcripts"""
    for cgp_id, tx in transcript_dict.iteritems():
        if 'jg' in tx.name2:
            continue  # this transcript was not assigned any parents
        ref_tx_ids = gene_tx_map[tx.name2]
        tx_seq = tx.get_cds(genome_fasta)
        for ref_tx_id in ref_tx_ids:
            if tx_biotype_map[ref_tx_id] != 'protein_coding':
                continue
            ref_tx = ref_transcript_dict[ref_tx_id]
            ref_tx_seq = ref_tx.get_cds(ref_genome_fasta)
            yield cgp_id, tx_seq, ref_tx_id, ref_tx_seq


def run_alignment_chunk(job, i, args, chunk):
    """
    Runs an alignment chunk through the blat/simpleChain pipeline
    :param i: # of chunk, for logging
    :param args: arguments to the pipeline
    :param chunk: List of (tx_id, tx_seq, ref_tx_id, ref_tx_seq) tuples
    :return: List of PSL output
    """
    job.fileStore.logToMaster('Beginning transcript alignment chunk {} for genome {}'.format(i, args['genome']))
    # temporary file paths
    ref_tmp_fasta = tools.fileOps.get_tmp_toil_file()
    tgt_tmp_fasta = tools.fileOps.get_tmp_toil_file()
    tmp_psl = tools.fileOps.get_tmp_toil_file()  # simpleChain does not like pipes
    blat_cmd = ['blat', '-noHead', '-extendThroughN', '-mask=lower', '-minIdentity=50', '-oneOff=1',
                ref_tmp_fasta, tgt_tmp_fasta, tmp_psl]
    chain_cmd = ['simpleChain', '-outPsl', tmp_psl, '/dev/stdout']
    results = []
    for tx_id, tx_seq, ref_tx_id, ref_tx_seq in chunk:
        tools.bio.write_fasta(ref_tmp_fasta, ref_tx_id, ref_tx_seq)
        tools.bio.write_fasta(tgt_tmp_fasta, tx_id, tx_seq)
        tools.procOps.run_proc(blat_cmd)
        chained_psl_strings = tools.procOps.call_proc_lines(chain_cmd)
        if len(chained_psl_strings) == 0:
            continue
        psl_recs = [tools.psl.PslRow(x.split()) for x in chained_psl_strings]
        best_psl_rec = sorted(psl_recs, key=lambda psl: psl.coverage)[-1]
        results.append(best_psl_rec.psl_string())
    return results


def merge(job, results, args):
    """
    Merge together chain files.
    :param results: list of promises from each alignment chunk
    :param args: arguments to the pipeline
    :return:
    """
    job.fileStore.logToMaster('Merging Alignment output for {}'.format(args['genome']), level=logging.INFO)
    tmp_results_file = tools.fileOps.get_tmp_file(tmp_dir=job.fileStore.getLocalTempDir())
    results_iter = itertools.chain.from_iterable(results)  # results is list of lists
    with open(tmp_results_file, 'w') as outf:
        for line in results_iter:
            outf.write(line + '\n')
    results_file_id = job.fileStore.writeGlobalFile(tmp_results_file)
    return results_file_id
