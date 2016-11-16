"""
Generate a hints database file from RNAseq alignments for AugustusTMR/AugustusCGP.

Expects a config file to be passed in with the paths to the files. Example:

[ANNOTATION]
annotation = /path/to/gff3

[INTRONBAM]
genome1 = /path/to/non_polyA_bam1.bam, /path/to/non_polyA_bam2.bam

[BAM]
genome1 = /path/to/fofn

The annotation field is optional, but will help AugustusCGP make better predictions.

BAM annotations can be put either under INTRONBAM or BAM. Any INTRONBAM will only have intron data loaded, and is
suitable for lower quality RNA-seq.

"""
import collections
import itertools
import logging
import os
import sqlite3
import argparse

import luigi
import luigi.contrib.sqla
import pyfasta
import pysam
from configobj import ConfigObj
from luigi.util import requires
from toil.common import Toil
from toil.job import Job

import tools.dataOps
import tools.fileOps
import tools.mathOps
import tools.misc
import tools.procOps
import tools.toilInterface
import tools.transcripts
import tools.hal
from base_tasks import HintsDbToilTask, HintsDbTask, HintsDbWrapperTask

logger = logging.getLogger(__name__)


class UserException(Exception):
    pass


class MissingFileException(UserException):
    pass


class BuildHints(HintsDbWrapperTask):
    """
    Main entry point. Parses input files, and launches the next task.
    """
    def parse_cfg(self):
        # configspec validates the input config file
        configspec = ['[ANNOTATION]', '__many__ = string', '[INTRONBAM]', '__many__ = list', '[BAM]', '__many__ = list']
        parser = ConfigObj(self.config, configspec=configspec)

        # convert the config into a new dict, parsing the BAMs
        cfg = collections.defaultdict(dict)
        target_genomes = set()
        if 'ANNOTATION' not in parser:
            cfg['ANNOTATION'] = {}
        else:
            for genome, annot in parser['ANNOTATION'].iteritems():
                annot = os.path.abspath(annot)
                if not os.path.exists(annot):
                    raise MissingFileException('Missing annotation file {}.'.format(annot))
                cfg['ANNOTATION'][genome] = annot
                target_genomes.add(genome)

        # if a given genome only has one BAM, it is a string. Fix this. Extract all paths from fofn files.
        for dtype in ['BAM', 'INTRONBAM']:
            if dtype not in parser:  # the user does not have to specify all field types
                cfg[dtype] = {}
                continue
            for genome in parser[dtype]:
                target_genomes.add(genome)
                path = parser[dtype][genome]
                if isinstance(path, str):
                    if not tools.misc.is_bam(path):
                        # this is a fofn
                        cfg[dtype][genome] = [os.path.abspath(x.rstrip()) for x in open(path)]
                    else:
                        # this is a single BAM
                        cfg[dtype][genome] = [os.path.abspath(path)]
                else:
                    cfg[dtype][genome] = [os.path.abspath(x) for x in path]
        # do some input validation
        for dtype in ['BAM', 'INTRONBAM']:
            for genome in cfg[dtype]:
                for bam in cfg[dtype][genome]:
                    if not os.path.exists(bam):
                        raise MissingFileException('Missing BAM {}.'.format(bam))
                    if not os.path.exists(bam + '.bai'):
                        raise MissingFileException('Missing BAM index {}.'.format(bam + '.bai'))
        for genome, annot in cfg['ANNOTATION'].iteritems():
            if not os.path.exists(annot):
                raise MissingFileException('Missing annotation file {}.'.format(annot))
        return cfg, tuple(target_genomes)

    def requires(self):
        cfg, target_genomes = self.parse_cfg()
        hint_paths = {}
        flat_fasta_paths = {}
        hal = os.path.abspath(self.hal)
        genomes = tools.hal.extract_genomes(hal)
        for genome in genomes:
            flat_fasta = self.clone(GenomeFlatFasta, genome=genome, cfg=cfg, hal=hal)
            yield flat_fasta
            flat_fasta_paths[genome] = flat_fasta.output().path
            if genome in target_genomes:
                annotation = cfg['ANNOTATION'].get(genome, None)
                hints = self.clone(GenerateHints, genome=genome, flat_fasta=flat_fasta.output().path,
                                   annotation=annotation, cfg=cfg)
                yield hints
                hint_paths[genome] = hints.output().path
        yield self.clone(BuildDb, cfg=cfg, hint_paths=hint_paths, flat_fasta_paths=flat_fasta_paths, genomes=genomes,
                         target_genomes=target_genomes)


class GenomeFlatFasta(HintsDbTask):
    """
    Flattens a genome fasta using pyfasta, copying it to the work directory. Requires the pyfasta package.
    """
    genome = luigi.Parameter()
    cfg = luigi.Parameter()
    hal = luigi.Parameter()

    def output(self):
        path = os.path.abspath(os.path.join(self.work_dir, self.genome + '.fa'))
        tools.fileOps.ensure_file_dir(path)
        return luigi.LocalTarget(path)

    def run(self):
        logger.info('Extracting fasta for {} from hal.'.format(self.genome))
        with self.output().open('w') as outf:
            cmd = ['hal2fasta', self.hal, self.genome]
            tools.procOps.run_proc(cmd, stdout=outf)
        logger.info('Flattening fasta for {}.'.format(self.genome))
        cmd = ['pyfasta', 'flatten', self.output().path]
        tools.procOps.run_proc(cmd)


@requires(GenomeFlatFasta)
class GenerateHints(HintsDbToilTask):
    """
    Generate hints for each genome as a separate Toil pipeline.
    """
    genome = luigi.Parameter()
    flat_fasta = luigi.Parameter()
    annotation = luigi.Parameter()
    cfg = luigi.Parameter()

    def output(self):
        hints = os.path.abspath(os.path.join(self.work_dir, self.genome + '.extrinsic_hints.gff'))
        return luigi.LocalTarget(hints)

    def run(self):
        logger.info('Beginning GenerateHints Toil pipeline for {}.'.format(self.genome))
        work_dir = os.path.abspath(os.path.join(self.work_dir, 'toil', self.genome))
        toil_options = self.prepare_toil_options(work_dir)
        toil_options.defaultMemory = '8G'  # TODO: don't hardcode this.
        generate_hints(self.genome, self.flat_fasta, self.cfg, self.annotation, self.output().path, toil_options)
        logger.info('Finished GenerateHints Toil pipeline for {}.'.format(self.genome))


class BuildDb(HintsDbTask):
    """
    Constructs the hints database from a series of reduced hints GFFs

    TODO: output() should be way smarter than this. Currently, it only checks if the indices have been created.
    """
    cfg = luigi.Parameter()
    hint_paths = luigi.Parameter()
    flat_fasta_paths = luigi.Parameter()
    genomes = luigi.TupleParameter()
    target_genomes = luigi.TupleParameter()

    def requires(self):
        for genome in self.target_genomes:
            flat_fasta = self.clone(GenomeFlatFasta, genome=genome, cfg=self.cfg)
            annotation = self.cfg['ANNOTATION'].get(genome, None)
            hints = self.clone(GenerateHints, genome=genome, flat_fasta=flat_fasta.output().path,
                               annotation=annotation, cfg=self.cfg)
            yield hints

    def output(self):
        tools.fileOps.ensure_file_dir(self.augustus_hints_db)
        return IndexTarget(self.augustus_hints_db)

    def run(self):
        for genome in self.genomes:
            logger.info('Loading sequence fpr {} into database.'.format(genome))
            base_cmd = ['load2sqlitedb', '--noIdx', '--clean', '--species={}'.format(genome),
                        '--dbaccess={}'.format(self.augustus_hints_db)]
            tools.procOps.run_proc(base_cmd + [self.flat_fasta_paths[genome]])
            if genome in self.hint_paths:
                logger.info('Loading hints for {} into database.'.format(genome))
                tools.procOps.run_proc(base_cmd + [self.hint_paths[genome]])
        logger.info('Indexing database.')
        cmd = ['load2sqlitedb', '--makeIdx', '--clean', '--dbaccess={}'.format(self.augustus_hints_db)]
        tools.procOps.run_proc(cmd)


class IndexTarget(luigi.Target):
    """
    luigi target that determines if the indices have been built on a hints database.
    """
    def __init__(self, db):
        self.db = db

    def exists(self, timeout=6000):
        con = sqlite3.connect(self.db, timeout=timeout)
        cur = con.cursor()
        r = []
        for idx in ['gidx', 'hidx']:
            query = 'PRAGMA index_info("{}")'.format(idx)
            try:
                v = cur.execute(query).fetchall()
            except sqlite3.OperationalError, exc:
                raise RuntimeError("query failed: {}\nOriginal error message: {}".format(query, exc))
            if len(v) > 0:
                r.append(v)
        return len(r) == 2


###
# Toil pipeline
###


def generate_hints(genome, flat_fasta, cfg, annotation, out_gff_path, toil_options):
    """
    Entry point for hints database Toil pipeline.
    """
    with Toil(toil_options) as toil:
        if not toil.options.restart:
            bam_file_ids = {'BAM': {}, 'INTRONBAM': {}}
            for dtype in ['BAM', 'INTRONBAM']:
                logger.info('Validating {}s for {}.'.format(dtype, genome))
                # validate BAMs
                fasta = pyfasta.Fasta(flat_fasta)
                fasta_sequences = {(x.split()[0], len(fasta[x])) for x in fasta.keys()}
                if genome not in cfg[dtype]:
                    continue
                for bam_path in cfg[dtype][genome]:
                    validate_bam_fasta_pairs(bam_path, fasta_sequences, genome)
                    is_paired = bam_is_paired(bam_path)
                    bam_file_ids[dtype][os.path.basename(bam_path)] = (toil.importFile('file://' + bam_path),
                                                                       toil.importFile('file://' + bam_path + '.bai'),
                                                                       is_paired)
                    is_paired_str = 'paired' if is_paired else 'not paired'
                    logger.info('BAM {} is valid and was inferred to be {}.'.format(os.path.basename(bam_path),
                                                                                    is_paired_str))
            input_file_ids = {'bams': bam_file_ids,
                              'annotation': toil.importFile('file://' + annotation) if annotation is not None else None}
            logger.info('{} has {} valid intron-only BAMs and {} valid BAMs. '
                        'Beginning Toil hints pipeline.'.format(genome, len(bam_file_ids['INTRONBAM']),
                                                                len(bam_file_ids['BAM'])))
            job = Job.wrapJobFn(setup_hints, input_file_ids)
            combined_hints = toil.start(job)
        else:
            logger.info('Restarting Toil hints pipeline for {}.'.format(genome))
            combined_hints = toil.restart()
        tools.fileOps.ensure_file_dir(out_gff_path)
        toil.exportFile(combined_hints, 'file://' + out_gff_path)


def setup_hints(job, input_file_ids):
    """
    Generates hints for a given genome with a list of BAMs. Will add annotation if it exists.

    Pipeline structure:
      filter_bam    cat_sort_bams    (build_intron_hints, build_exon_hints)
        ^  V           ^ V                 ^  V
    setup_hints -> merge_bams -------> build_hints -------> cat_hints
          V generate_annotation_hints ------------------------^

    Each main step (filter_bam, cat_sort_bams, build_intron_hints, build_exon_hints) are done on a subset of references
    that are then combined at the cat_hints step.
    """
    filtered_bam_file_ids = {'BAM': collections.defaultdict(list), 'INTRONBAM': collections.defaultdict(list)}
    for dtype, bam_dict in input_file_ids['bams'].iteritems():
        if len(bam_dict) == 0:
            continue
        # Since BAMs are valid, we can assume that they all share the same header
        bam_file_id, bai_file_id, is_paired = bam_dict.values()[0]
        sam_handle = pysam.Samfile(job.fileStore.readGlobalFile(bam_file_id))
        # generate reference grouping that will be used downstream until final cat step
        grouped_references = [tuple(x) for x in group_references(sam_handle)]
        for original_path, (bam_file_id, bai_file_id, is_paired) in bam_dict.iteritems():
            for reference_subset in grouped_references:
                j = job.addChildJobFn(filter_bam, bam_file_id, bai_file_id, reference_subset, is_paired)
                filtered_bam_file_ids[dtype][reference_subset].append(j.rv())
    if input_file_ids['annotation'] is not None:
        j = job.addChildJobFn(generate_annotation_hints, input_file_ids['annotation'])
        annotation_hints_file_id = j.rv()
    else:
        annotation_hints_file_id = None
    # returns path to filtered GFF output
    return job.addFollowOnJobFn(merge_bams, filtered_bam_file_ids, annotation_hints_file_id).rv()


def filter_bam(job, bam_file_id, bai_file_id, reference_subset, is_paired):
    """
    Slices out a chromosome from a BAM, re-sorts by name, filters the reads, then re-sorts by position.
    filterBam does weird things when piped to stdout, so I don't do that.
    """
    tmp_filtered = tools.fileOps.get_tmp_toil_file(suffix='filtered.bam')
    bam_path = job.fileStore.readGlobalFile(bam_file_id)
    job.fileStore.readGlobalFile(bai_file_id, bam_path + '.bai')
    sort_tmp = tools.fileOps.get_tmp_toil_file()
    cmd = [['samtools', 'view', '-b', bam_path],
           ['samtools', 'sort', '-O', 'bam', '-T', sort_tmp, '-n', '-l', '0', '-'],
           ['filterBam', '--uniq', '--in', '/dev/stdin', '--out', tmp_filtered]]
    if is_paired is True:
        cmd[-1].extend(['--paired', '--pairwiseAlignments'])
    cmd[0].extend(reference_subset)
    tools.procOps.run_proc(cmd)
    out_filter = tools.fileOps.get_tmp_toil_file(suffix='sorted.filtered.bam')
    sort_cmd = ['samtools', 'sort', '-O', 'bam', '-T', sort_tmp, tmp_filtered]
    tools.procOps.run_proc(sort_cmd, stdout=out_filter)
    filtered_bam_file_id = job.fileStore.writeGlobalFile(out_filter)
    return filtered_bam_file_id


def merge_bams(job, filtered_bam_file_ids, annotation_hints_file_id):
    """
    Takes a dictionary mapping reference chunks to filtered BAMs. For each reference chunk, these BAMs will be
    first concatenated then sorted, then passed off to hint building.
    """
    merged_bam_file_ids = {'BAM': {}, 'INTRONBAM': {}}
    for dtype in filtered_bam_file_ids:
        for ref_group, file_ids in filtered_bam_file_ids[dtype].iteritems():
            merged_bam_file_ids[dtype][ref_group] = job.addChildJobFn(cat_sort_bams, file_ids).rv()
    return job.addFollowOnJobFn(build_hints, merged_bam_file_ids, annotation_hints_file_id).rv()


def cat_sort_bams(job, bam_file_ids):
    """
    Takes a list of bam file IDs and combines/sorts them.

    TODO: the 4096 file hack below is hacky. Should only be a problem for very fragmented references.
    """
    bamfiles = [job.fileStore.readGlobalFile(x) for x in bam_file_ids]
    # cat only 4095 bams at a time to avoid bash command length problems
    catfile = tools.fileOps.get_tmp_toil_file()
    sam_iter = tools.dataOps.grouper(bamfiles, 4095)
    # do the first one
    cmd = ['samtools', 'cat', '-o', catfile]
    cmd.extend(sam_iter.next())
    tools.procOps.run_proc(cmd)
    # do any subsequent ones left, creating a new file each time
    for more in sam_iter:
        old_catfile = catfile
        catfile = tools.fileOps.get_tmp_toil_file()
        cmd = ['samtools', 'cat', '-o', catfile, old_catfile]
        cmd.extend(more)
        tools.procOps.run_proc(cmd)
    # combine and merge
    merged = tools.fileOps.get_tmp_toil_file()
    sort_tmp = tools.fileOps.get_tmp_toil_file()
    cmd = ['samtools', 'sort', '-O', 'bam', '-T', sort_tmp, catfile]
    tools.procOps.run_proc(cmd, stdout=merged)
    return job.fileStore.writeGlobalFile(merged)


def build_hints(job, merged_bam_file_ids, annotation_hints_file_id):
    """
    Takes the merged BAM for a genome and produces both intron and exon hints.
    """
    intron_hints_file_ids = []
    exon_hints_file_ids = []
    for dtype in merged_bam_file_ids:
        for ref_group, file_ids in merged_bam_file_ids[dtype].iteritems():
            intron_hints_file_ids.append(job.addChildJobFn(build_intron_hints, file_ids).rv())
            if dtype == 'BAM':
                exon_hints_file_ids.append(job.addChildJobFn(build_exon_hints, file_ids).rv())
    return job.addFollowOnJobFn(cat_hints, intron_hints_file_ids, exon_hints_file_ids, annotation_hints_file_id).rv()


def build_intron_hints(job, merged_bam_file_id):
    """Builds intronhints from a BAM. Returns a fileID to the hints."""
    bam_file = job.fileStore.readGlobalFile(merged_bam_file_id)
    intron_gff_path = tools.fileOps.get_tmp_toil_file()
    cmd = ['bam2hints', '--intronsonly', '--in', bam_file, '--out', intron_gff_path]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(intron_gff_path)


def build_exon_hints(job, merged_bam_file_id):
    """Builds exonhints from a BAM Returns a fileID to the hints."""
    bam_file = job.fileStore.readGlobalFile(merged_bam_file_id)
    cmd = [['bam2wig', bam_file],
           ['wig2hints.pl', '--width=10', '--margin=10', '--minthresh=2', '--minscore=4', '--prune=0.1', '--src=W',
            '--type=ep', '--UCSC=/dev/null', '--radius=4.5', '--pri=4', '--strand=.']]
    exon_gff_path = tools.fileOps.get_tmp_toil_file()
    tools.procOps.run_proc(cmd, stdout=exon_gff_path)
    return job.fileStore.writeGlobalFile(exon_gff_path)


def generate_annotation_hints(job, annotation_hints_file_id):
    """
    Converts the annotation file into hints. First converts the gff3 directly to genePred so we can make use
    of the transcript library.

    Hints are derived from both CDS exonic intervals and intron intervals
    """
    annotation_gff3 = job.fileStore.readGlobalFile(annotation_hints_file_id)
    tm_gp = tools.fileOps.get_tmp_toil_file()
    cmd = ['gff3ToGenePred', '-rnaNameAttr=transcript_id', '-geneNameAttr=gene_id', '-honorStartStopCodons',
           annotation_gff3, tm_gp]
    tools.procOps.run_proc(cmd)
    tx_dict = tools.transcripts.get_gene_pred_dict(tm_gp)
    hints = []
    for tx_id, tx in tx_dict.iteritems():
        if tx.cds_size == 0:
            continue
        # rather than try to re-do the arithmetic, we will use the get_bed() function to convert this transcript
        cds_tx = tools.transcripts.Transcript(tx.get_bed(new_start=tx.thick_start, new_stop=tx.thick_stop))
        for intron in cds_tx.intron_intervals:
            r = [intron.chromosome, 'a2h', 'intron', intron.start + 1, intron.stop, 0, intron.strand, '.',
                 'grp={};src=M;pri=2'.format(tx_id)]
            hints.append(r)
        for exon in cds_tx.exon_intervals:
            r = [exon.chromosome, 'a2h', 'CDS', exon.start + 1, exon.stop, 0, exon.strand, '.',
                 'grp={};src=M;pri=2'.format(tx_id)]
            hints.append(r)
    annotation_hints_gff = tools.fileOps.get_tmp_toil_file()
    tools.fileOps.print_rows(annotation_hints_gff, hints)
    return job.fileStore.writeGlobalFile(annotation_hints_gff)


def cat_hints(job, intron_hints_file_ids, exon_hints_file_ids, annotation_hints_file_id):
    """Returns file ID to combined, sorted hints"""
    cat_hints = tools.fileOps.get_tmp_toil_file()
    with open(cat_hints, 'w') as outf:
        for file_id in itertools.chain(intron_hints_file_ids, exon_hints_file_ids):
            f = job.fileStore.readGlobalFile(file_id)
            for line in open(f):
                outf.write(line)
        if annotation_hints_file_id is not None:
            f = job.fileStore.readGlobalFile(annotation_hints_file_id)
            for line in open(f):
                outf.write(line)
    cmd = [['sort', '-n', '-k4,4', cat_hints],
           ['sort', '-s', '-n', '-k5,5'],
           ['sort', '-s', '-n', '-k3,3'],
           ['sort', '-s', '-k1,1'],
           ['join_mult_hints.pl']]
    combined_hints = tools.fileOps.get_tmp_toil_file()
    tools.procOps.run_proc(cmd, stdout=combined_hints)
    return job.fileStore.writeGlobalFile(combined_hints)


###
# Functions
###


def validate_bam_fasta_pairs(bam_path, fasta_sequences, genome):
    """
    Make sure that this BAM is actually aligned to this fasta. Every sequence should be the same length. Sequences
    can exist in the reference that do not exist in the BAM, but not the other way around.
    """
    handle = pysam.Samfile(bam_path, 'rb')
    bam_sequences = {(n, s) for n, s in zip(*[handle.references, handle.lengths])}
    difference = bam_sequences - fasta_sequences
    if len(difference) > 0:
        base_err = 'Error: BAM {} has the following sequence/length pairs not found in the {} fasta: {}.'
        err = base_err.format(bam_path, genome, ','.join(['-'.join(map(str, x)) for x in difference]))
        raise UserException(err)
    missing_seqs = fasta_sequences - bam_sequences
    if len(missing_seqs) > 0:
        base_msg = 'BAM {} does not have the following sequence/length pairs in its header: {}'
        msg= base_msg.format(bam_path, ','.join(['-'.join(map(str, x)) for x in missing_seqs]))
        logger.warning(msg)


def bam_is_paired(bam_path, num_reads=20000, paired_cutoff=0.75):
    """
    Infers the paired-ness of a bam file.
    """
    sam = pysam.Samfile(bam_path)
    count = 0
    for rec in itertools.islice(sam, num_reads):
        if rec.is_paired:
            count += 1
    if tools.mathOps.format_ratio(count, num_reads) > 0.75:
        return True
    elif tools.mathOps.format_ratio(count, num_reads) < 1 - paired_cutoff:
        return False
    else:
        raise UserException("Unable to infer pairing from bamfile {}".format(bam_path))


def group_references(sam_handle, num_bases=10 ** 7, max_seqs=1000):
    """
    Group up references by num_bases, unless that exceeds max_seqs. A greedy implementation of the bin packing problem.
    """
    name_iter = itertools.izip(*[sam_handle.references, sam_handle.lengths])
    name, size = name_iter.next()
    this_bin = [name]
    bin_base_count = size
    num_seqs = 1
    for name, size in name_iter:
        bin_base_count += size
        num_seqs += 1
        if bin_base_count >= num_bases or num_seqs > max_seqs:
            yield this_bin
            this_bin = [name]
            bin_base_count = size
            num_seqs = 1
        else:
            this_bin.append(name)
    yield this_bin


###
# Entry point without using luigi
###


def parse_args():
    """
    Provides the ability to run directly from this script, bypassing the luigi wrapper. If you go this route, you
    cannot control the number of concurrent toil pipelines.
    :return: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--hal', required=True)
    parser.add_argument('--augustus-hints-db', default='augustus_hints.db')
    parser.add_argument('--work-dir', default='./hints_work')
    # parallelism
    parser.add_argument('--workers', default=10)
    # toil options
    parser.add_argument('--batchSystem', default='singleMachine')
    parser.add_argument('--maxCores', default=16, type=int)
    parser.add_argument('--logLevel', default='WARNING')
    parser.add_argument('--cleanWorkDir', default='onSuccess')
    parser.add_argument('--parasolCommand', default=None)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    workers = args.workers
    del args.workers  # hack because workers is actually a argument to luigi, not RunCat
    luigi.build([BuildHints(**vars(args))], logging_conf_file='logging.cfg', workers=workers)
