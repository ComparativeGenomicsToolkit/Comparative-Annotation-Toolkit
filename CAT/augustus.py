"""
Runs AugustusTM(R) on the input transMap records.

This program takes as input a genePred of transMap output, an attributes table describing the biotype and
gene-transcript relationships, and an optional RNAseq hints database and runs Augustus on each transcript.

Transcripts are evaluated for the following features to help decide how they are converted into hints:

    1) OriginalIntrons. If any new gaps are not within a wiggle distance (in transcript space) of original introns,
        do not provide them as hints to Augustus.
    2) Original thick_start/thick_stop. If we did not map over the original start/stop, do not provide them as hints.
    3) Original tss/tts. If we did not map over the original transcription start/stop, do not provide them as hints.

"""
import logging
import os

import sqlalchemy
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import sessionmaker
from toil.job import Job
from toil.common import Toil

import tools.bio
import tools.dataOps
import tools.fileOps
import tools.procOps
import tools.psl
import tools.tm2hints
import tools.transcripts


def reflect_hints_db(db_path):
    """
    Reflect the database schema of the hints database, automapping the existing tables
    :param db_path: path to hints sqlite database
    :return: sqlalchemy.MetaData object, sqlalchemy.orm.Session object
    """
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_path))
    metadata = sqlalchemy.MetaData()
    metadata.reflect(bind=engine)
    Base = automap_base(metadata=metadata)
    Base.prepare()
    speciesnames = Base.classes.speciesnames
    seqnames = Base.classes.seqnames
    hints = Base.classes.hints
    featuretypes = Base.classes.featuretypes
    Session = sessionmaker(bind=engine)
    session = Session()
    return speciesnames, seqnames, hints, featuretypes, session


def get_rnaseq_hints(genome, chromosome, start, stop, speciesnames, seqnames, hints, featuretypes, session):
    """
    Extracts RNAseq hints from RNAseq hints database.
    :param genome: genome (table) to query
    :param chromosome: Chromosome to extract information from
    :param start: start position on chromosome
    :param stop: stop position in chromosome
    :param speciesnames: speciesnames Table from reflect_hints_db
    :param seqnames: seqnames Table from reflect_hints_db
    :param hints: hints Table from reflect_hints_db
    :param featuretypes: featuretypes Table from reflect_hints_db
    :param session: Session object from reflect_hints_db
    :return: GFF formatted string.
    """
    speciesid = session.query(speciesnames.speciesid).filter_by(speciesname=genome)
    seqnr = session.query(seqnames.seqnr).filter(
        sqlalchemy.and_(
            seqnames.speciesid.in_(speciesid),
            (seqnames.seqname == chromosome)))
    query = session.query(hints, featuretypes).filter(
            sqlalchemy.and_(
                hints.speciesid.in_(speciesid),
                hints.seqnr.in_(seqnr),
                hints.start >= start,
                hints.end <= stop,
                featuretypes.typeid == hints.type))
    hints = []
    for h, f in query:
        tags = 'pri=3;src={};mult={}'.format(h.esource, h.mult)
        l = [chromosome, h.source, f.typename, h.start, h.end, h.score, '.', '.', tags]
        hints.append('\t'.join(map(str, l)) + '\n')
    return '\n'.join(hints)


def rename_transcripts(aug_output, valid_txs, cfg_version, name):
    """
    Renames overlapping transcripts augIX-ID, where X is the index of the extrinsic.cfg file, e.g. 1 or 2 and where
    ID is the transmap alignment ID, use augIX-n-ID if more than 1 transcript overlaps the alignment
    formats these transcripts into a string
    """
    name_map = {}
    for i, x in enumerate(valid_txs):
        if i > 0:
            name_map[x] = "augI{}-{}-{}".format(cfg_version, i + 1, name)
        else:
            name_map[x] = "augI{}-{}".format(cfg_version, name)
    print name_map
    features = {"exon", "CDS", "start_codon", "stop_codon", "tts", "tss"}
    for x in aug_output:
        if x.startswith("#"):
            continue
        if "AUGUSTUS" in x:
            x = x.split("\t")
            if x[2] in features:
                t = x[-1].split()
                n = t[-3].split('"')[1]
                if n not in name_map:
                    continue  # skip transcripts filtered out previously
                try:
                    t[-1] = t[-3] = '"{}";'.format(name_map[n])
                except:
                    assert False, (name_map, aug_output, name)
                t = " ".join(t)
                x[-1] = t
                yield x


def setup(job, args, coding_gp, chunk_size=5):
    """
    Entry function for running AugustusTM(R). Loads the genome fasta into the fileStore then spins up chunks of
    jobs.
    :param args: args from Luigi pipeline
    :param coding_gp: genePred with only coding transcripts
    :param chunk_size: number of transcripts to process per-job
    :return:
    """
    job.fileStore.logToMaster('Beginning Augustus run on {}'.format(args['genome']), level=logging.INFO)
    fasta_file_id = job.fileStore.writeGlobalFile(args['genome_fasta'])
    # assume that this fasta has been flattened
    gdx_file_id = job.fileStore.writeGlobalFile(args['genome_fasta'] + '.gdx')
    flat_file_id = job.fileStore.writeGlobalFile(args['genome_fasta'] + '.flat')
    tm_cfg_file_id = job.fileStore.writeGlobalFile(os.path.abspath(args['tm_cfg']))
    if args['augustus_hints_db'] is not None:
        hints_db_file_id = job.fileStore.writeGlobalFile(args['augustus_hints_db'])
        tmr_cfg_file_id = job.fileStore.writeGlobalFile(args['tmr_cfg'])
    else:
        hints_db_file_id = tmr_cfg_file_id = None
    # load the reference and transMap PSL data necessary for filtering transMap hints
    ref_psl_dict = tools.psl.get_alignment_dict(args['ref_psl'])
    tm_psl_dict = tools.psl.get_alignment_dict(args['tm_psl'])
    job.fileStore.logToMaster('Augustus loaded reference and transMap PSLs')
    results = []
    for i, chunk in enumerate(tools.dataOps.grouper(tools.transcripts.gene_pred_iterator(coding_gp), chunk_size)):
        if i > 40:
            break
        grouped_recs = {}
        for tx_id, tx in chunk:
            grouped_recs[tx_id] = [tx, tm_psl_dict[tx_id], ref_psl_dict[tools.psl.remove_alignment_number(tx_id)]]
        j = job.addChildJobFn(run_augustus_chunk, i, args, grouped_recs, fasta_file_id, gdx_file_id, tm_cfg_file_id,
                              flat_file_id, hints_db_file_id, tmr_cfg_file_id)
        results.append(j.rv())
    return results


def run_augustus_chunk(job, i, args, grouped_recs, fasta_file_id, gdx_file_id, tm_cfg_file_id, flat_file_id,
                       hints_db_file_id, tmr_cfg_file_id, padding=20000):
    """
    Runs augustus on a chunk of genePred objects.
    :param i: chunk ID. for logging.
    :param args: Arguments passed by Luigi
    :param grouped_recs: Chunk of (tx_id, GenePredTranscript) tuples
    :param fasta_file_id: fileStore ID for the genome fasta
    :param tm_cfg_file_id: fileStore ID for the AugustusTM config
    :param flat_file_id: fileStore ID for the genome fasta flat file
    :param hints_db_file_id: fileStore ID for the hintsDb, if it exists
    :param tmr_cfg_file_id: fileStore ID for the AugustusTMR config, if the hints db exists
    :param padding: Amount of padding to put around the transcript region, allowing Augustus to add 5'/3' sequence.
    :return: Augustus output for this chunk
    """
    job.fileStore.logToMaster('Beginning chunk {}'.format(i))
    # we have to explicitly place fasta, flat file and  gdx with the correct naming scheme for pyfasta
    work_dir = job.fileStore.getLocalTempDir()
    fasta_local_path = os.path.join(work_dir, 'genome.fasta')
    gdx_local_path = os.path.join(work_dir, 'genome.fasta.gdx')
    flat_local_path = os.path.join(work_dir, 'genome.fasta.flat')
    job.fileStore.readGlobalFile(fasta_file_id, fasta_local_path)
    job.fileStore.readGlobalFile(gdx_file_id, gdx_local_path)
    job.fileStore.readGlobalFile(flat_file_id, flat_local_path)
    fasta = tools.bio.get_sequence_dict(fasta_local_path, upper=False)  # maintain softmasking
    job.fileStore.logToMaster('Chunk {} successfully loaded the fasta'.format(i))
    tm_cfg_file = job.fileStore.readGlobalFile(tm_cfg_file_id)

    if hints_db_file_id is not None:  # we are running TMR mode as well
        job.fileStore.logToMaster('Chunk {} is in TMR mode'.format(i))
        tmr_cfg_file = job.fileStore.readGlobalFile(tmr_cfg_file_id)
        hints_db_file = job.fileStore.readGlobalFile(hints_db_file_id)
        speciesnames, seqnames, hints, featuretypes, session = reflect_hints_db(hints_db_file)
        job.fileStore.logToMaster('Chunk {} successfully loaded the hints database'.format(i))
    else:
        job.fileStore.logToMaster('Chunk {} is in TM mode'.format(i))

    # start iteratively running Augustus on this chunk
    results = []
    for tm_tx, tm_psl, ref_psl in grouped_recs.itervalues():
        if len(tm_tx) > 3 * 10 ** 6:  # no huge transcripts
            continue
        tm_hints = tools.tm2hints.tm_to_hints(tm_tx, ref_psl, tm_psl)
        chromosome = tm_tx.chromosome
        start = max(tm_tx.start - padding, 0)
        stop = min(tm_tx.stop + padding, len(fasta[chromosome]))
        if hints_db_file_id is not None:
            rnaseq_hints = get_rnaseq_hints(args['genome'], chromosome, start, stop, speciesnames, seqnames, hints,
                                            featuretypes, session)
            hint = ''.join([tm_hints, rnaseq_hints])
            results.append(run_augustus(job, hint, fasta, tm_tx, tmr_cfg_file, start, stop, cfg_version=2))
        else:
            hint = tm_hints
        transcripts = run_augustus(job, hint, fasta, tm_tx, tm_cfg_file, start, stop, cfg_version=1)
        if len(transcripts) > 0:  # we may not have found anything
            results.append(transcripts)
    return results


def run_augustus(job, hint, fasta, tm_tx, cfg_file, start, stop, cfg_version=1):
    """
    Runs Augustus.
    :param job: job instance
    :param hint: GFF formatted hint string
    :param fasta: Pyfasta object
    :param tm_tx: GenePredTranscript object
    :param cfg_file: config file
    :param cfg_version: config file version
    :return: GFF formatted output from Augustus
    """
    with tools.fileOps.TemporaryFilePath(tmp_dir=job.fileStore.getLocalTempDir()) as hints_out, \
            tools.fileOps.TemporaryFilePath(tmp_dir=job.fileStore.getLocalTempDir()) as fasta_out:
        with open(hints_out, 'w') as hints_out_handle, open(fasta_out, 'w') as fasta_out_handle:
            hints_out_handle.write(hint)
            fasta_out_handle.write(fasta[tm_tx.chromosome][start:stop])
            cmd = ['augustus', fasta_out, '--predictionStart=-{}'.format(start), '--predictionEnd=-{}'.format(stop),
                   '--extrinsicCfgFile={}'.format(cfg_file), '--hintsfile={}'.format(hints_out), '--UTR=on',
                   '--alternatives-from-evidence=0', '--species=human', '--allow_hinted_splicesites=atac',
                   '--protein=0', '--softmasking=1']
        aug_output = tools.procOps.call_proc(cmd)
    aug_output = aug_output.split("\n")
    print 'aug output:\n{}\ntx:{}'.format(aug_output, tm_tx.name)
    # extract the transcript lines
    tx_entries = [x.split() for x in aug_output if "\ttranscript\t" in x]
    # filter out transcripts that do not overlap the alignment range
    valid_txs = [x[-1] for x in tx_entries if not (int(x[4]) < tm_tx.start or int(x[3]) > tm_tx.stop)]
    # rename transcript based on cfg version, and make names unique
    transcripts = list(rename_transcripts(aug_output, valid_txs, cfg_version, tm_tx.name))
    print 'transcripts after renaming:\n{}'.format(transcripts)
    return transcripts


def augustus(args, coding_gp, toil_options):
    toil_options.cleanWorkDir = 'never'
    j = Job.wrapJobFn(setup, args, coding_gp)
    results = Job.Runner.startToil(j, toil_options)
    return results
