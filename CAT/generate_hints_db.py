"""
Generate a hints database file from RNAseq alignments for AugustusTMR/AugustusCGP.

Expects a config file to be passed in with the paths to the files. Example:

[FASTA]
genome1 = /path/to/fasta

[ANNOTATION]
annotation = /path/to/gff3

[BAM]
genome1 = /path/to/bam1.bam, /path/to/bam2.bam

The annotation field is optional, but will help AugustusCGP make better predictions.

"""
import collections
import itertools
import logging
import os
import shutil
import sqlite3
import tempfile

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

logger = logging.getLogger(__name__)


class UserException(Exception):
    pass


class HintsDbTask(luigi.Task):
    # path to config file
    config = luigi.Parameter(default='hints_cfg.cfg')
    augustus_hints_db = luigi.Parameter(default='augustus_hints.db')
    work_dir = luigi.Parameter(default=os.path.join(tempfile.gettempdir(), __name__))
    # Toil options
    batchSystem = luigi.Parameter(default='singleMachine', significant=False)
    maxCores = luigi.IntParameter(default=16, significant=False)
    logLevel = luigi.Parameter(default='WARNING', significant=False)  # this is passed to toil
    cleanWorkDir = luigi.Parameter(default='onSuccess', significant=False)  # debugging option
    parasolCommand = luigi.Parameter(default=None, significant=False)
    defaultMemory = luigi.IntParameter(default=8 * 1024 ** 3, significant=False)

    def __repr__(self):
        """override the repr to make logging cleaner"""
        return 'HintsDbTask: {}'.format(self.__class__.__name__)


class HintsDbWrapperTask(HintsDbTask, luigi.WrapperTask):
    """add WrapperTask functionality to PipelineTask"""
    pass


class HintsDbToilTask(HintsDbTask):
    """
    Task for launching toil pipelines from within luigi.
    """
    resources = {'toil': 1}  # all toil pipelines use 1 toil

    def prepare_toil_options(self, work_dir):
        """
        Prepares a Namespace object for Toil which has all defaults, overridden as specified
        Will see if the jobStore path exists, and if it does, assume that we need to add --restart
        :param work_dir: Parent directory where toil work will be done. jobStore will be placed inside. Will be used
        to fill in the workDir class variable.
        :return: Namespace
        """
        job_store = os.path.join(work_dir, 'jobStore')
        tools.fileOps.ensure_file_dir(job_store)
        toil_args = self.get_toil_defaults()
        toil_args.__dict__.update(vars(self))
        if os.path.exists(job_store):
            try:
                root_job = open(os.path.join(job_store, 'rootJobStoreID')).next().rstrip()
                if not os.path.exists(os.path.join(job_store, 'tmp', root_job)):
                    shutil.rmtree(job_store)
                else:
                    toil_args.restart = True
            except OSError:
                toil_args.restart = True
            except IOError:
                shutil.rmtree(job_store)
        job_store = 'file:' + job_store
        toil_args.jobStore = job_store
        return toil_args

    def get_toil_defaults(self):
        """
        Extracts the default toil options as a dictionary, setting jobStore to None
        :return: dict
        """
        parser = Job.Runner.getDefaultArgumentParser()
        namespace = parser.parse_args([''])  # empty jobStore attribute
        namespace.jobStore = None  # jobStore attribute will be updated per-batch
        return namespace

    def __repr__(self):
        """override the PipelineTask repr to report the batch system being used"""
        base_repr = super(HintsDbToilTask, self).__repr__()
        return 'Toil' + base_repr + ' using batchSystem {}'.format(self.batchSystem)


class BuildHints(HintsDbWrapperTask):
    """
    Main entry point. Parses input files, and launches the next task.
    """
    def parse_cfg(self):
        # configspec validates the input config file
        configspec = ['[FASTA]', '__many__ = string', '[ANNOTATION]', '__many__ = string', '[BAM]', '__many__ = list']
        parser = ConfigObj(self.config, configspec=configspec)
        # if a given genome only has one BAM, it is a string. Fix this.
        for genome in parser['BAM']:
            path = parser['BAM'][genome]
            if isinstance(path, str):
                if not tools.misc.is_bam(path):
                    # this is a fofn
                    parser['BAM'][genome] = [x.rstrip() for x in open(path)]
                else:
                    # this is a single BAM
                    parser['BAM'][genome] = [path]
        # do some input validation
        for genome in parser['BAM']:
            for bam in parser['BAM'][genome]:
                assert os.path.exists(bam)
                assert os.path.exists(bam + '.bai')
        # TODO: validate all inputs. Some validation happens in the toil pipeline as well.
        return parser

    def requires(self):
        cfg = self.parse_cfg()
        hint_paths = {}
        flat_fasta_paths = {}
        for genome in cfg['FASTA']:
            flat_fasta = self.clone(GenomeFlatFasta, genome=genome, cfg=cfg)
            yield flat_fasta
            flat_fasta_paths[genome] = flat_fasta.output().path
            annotation = cfg['ANNOTATION'].get(genome, None)
            hints = self.clone(GenerateHints, genome=genome, flat_fasta=flat_fasta.output().path,
                               annotation=annotation, cfg=cfg)
            yield hints
            hint_paths[genome] = hints.output().path
        yield self.clone(BuildDb, cfg=cfg, hint_paths=hint_paths, flat_fasta_paths=flat_fasta_paths)


class GenomeFlatFasta(HintsDbTask):
    """
    Flattens a genome fasta using pyfasta, copying it to the work directory. Requires the pyfasta package.
    """
    genome = luigi.Parameter()
    cfg = luigi.Parameter()

    def output(self):
        path = os.path.abspath(os.path.join(self.work_dir, self.genome + '.fasta'))
        tools.fileOps.ensure_file_dir(path)
        return luigi.LocalTarget(path)

    def run(self):
        logger.info('Flattening fasta for {}.'.format(self.genome))
        shutil.copyfile(self.cfg['FASTA'][self.genome], self.output().path)
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
    resources = {'toil': 1}  # all toil pipelines use 1 toil

    def output(self):
        hints = os.path.abspath(os.path.join(self.work_dir, self.genome + '.reduced_hints.gff'))
        return luigi.LocalTarget(hints)

    def run(self):
        logger.info('Beginning GenerateHints toil pipeline for {}.'.format(self.genome))
        work_dir = os.path.abspath(os.path.join(self.work_dir, 'toil', self.genome))
        toil_options = self.prepare_toil_options(work_dir)
        toil_options.defaultMemory = '8G'  # TODO: don't hardcode this.
        completed = generate_hints(self.genome, self.flat_fasta, self.cfg, self.annotation, self.output().path,
                                   toil_options)
        if completed is False:  # we did not have any hints to generate for this genome
            logger.info('Did not generate hints for {} due to a lack of BAMs/annotation'.format(self.genome))
            self.output().open('w').close()
        else:
            logger.info('Finished generating hints for {}.'.format(self.genome))


class BuildDb(HintsDbTask):
    """
    Constructs the hints database from a series of reduced hints GFFs

    TODO: output() should be way smarter than this. Currently, it only checks if the indices have been created.
    """
    cfg = luigi.Parameter()
    hint_paths = luigi.Parameter()
    flat_fasta_paths = luigi.Parameter()

    def requires(self):
        for genome in self.cfg['FASTA']:
            flat_fasta = self.clone(GenomeFlatFasta, genome=genome, cfg=self.cfg)
            annotation = self.cfg['ANNOTATION'].get(genome, None)
            hints = self.clone(GenerateHints, genome=genome, flat_fasta=flat_fasta.output().path,
                               annotation=annotation, cfg=self.cfg)
            yield hints

    def output(self):
        tools.fileOps.ensure_file_dir(self.augustus_hints_db)
        return IndexTarget(self.augustus_hints_db)

    def run(self):
        for genome in self.cfg['FASTA']:
            logger.info('Loading finished hints for {} into database.'.format(genome))
            base_cmd = ['load2sqlitedb', '--noIdx', '--clean', '--species={}'.format(genome),
                        '--dbaccess={}'.format(self.augustus_hints_db)]
            tools.procOps.run_proc(base_cmd + [self.flat_fasta_paths[genome]])
            if genome in self.hint_paths and os.path.getsize(self.hint_paths[genome]) > 0:
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
    if genome in cfg['BAM'] and cfg['BAM'][genome] is not None:
        logger.info('Validating BAMs for {}.'.format(genome))
        # validate BAMs
        fasta = pyfasta.Fasta(flat_fasta)
        fasta_sequences = {x.split()[0] for x in fasta.keys()}
        for bam_path in cfg['BAM'][genome]:
            validate_bam_fasta_pairs(bam_path, fasta_sequences, genome)
            logger.info('BAM {} for {} is valid.'.format(bam_path, genome))
        logger.info('All BAMs valid for {}, beginning Toil hints pipeline.'.format(genome))
    # start pipeline
    with Toil(toil_options) as toil:
        if not toil.options.restart:
            if genome in cfg['BAM'] and cfg['BAM'][genome] is not None:
                bam_file_ids = {}
                for bam_path in cfg['BAM'][genome]:
                    is_paired = bam_is_paired(bam_path)
                    bam_file_ids[os.path.basename(bam_path)] = (toil.importFile('file://' + bam_path),
                                                                toil.importFile('file://' + bam_path + '.bai'),
                                                                is_paired)
                    is_paired_str = 'paired' if is_paired else 'not paired'
                    logger.info('BAM {} was inferred to be {}.'.format(os.path.basename(bam_path), is_paired_str))
            else:
                bam_file_ids = None
            input_file_ids = {'bams': bam_file_ids,
                              'annotation': toil.importFile('file://' + annotation) if annotation is not None else None}
            job = Job.wrapJobFn(setup_hints, input_file_ids, genome)
            combined_hints = toil.start(job)
        else:
            combined_hints = toil.restart()
        if combined_hints is not None:
            tools.fileOps.ensure_file_dir(out_gff_path)
            toil.exportFile(combined_hints, 'file://' + out_gff_path)
            return True
        else:
            return False


def setup_hints(job, input_file_ids, genome):
    """
    Generates hints for a given genome with a list of BAMs. Will add annotation if it exists.

    Pipeline structure:
      filter_bam    cat_sort_bams    (build_intron_hints, build_exon_hints)
        ^  V           ^ V                 ^  V
    setup_hints -> merge_bams -------> build_hints -------> cat_hints
          V generate_annotation_hints ------------------------^

    Each main step (filter_bam, cat_sort_bams, build_intron_hints, build_exon_hints) are done on a subset of references
    and then combined at the cat step.
    """
    job.fileStore.logToMaster('Beginning hints production for {}'.format(genome), level=logging.INFO)
    # If we have no BAMs, we skip this step
    if input_file_ids['bams'] is not None:
        # Since BAMs are valid, we can assume that they all share the same header
        bam_file_id, bai_file_id, is_paired = input_file_ids['bams'].values()[0]
        sam_handle = pysam.Samfile(job.fileStore.readGlobalFile(bam_file_id))
        # generate reference grouping that will be used downstream until final cat step
        grouped_references = [tuple(x) for x in group_references(sam_handle)]
        filtered_bam_file_ids = collections.defaultdict(list)
        for original_path, (bam_file_id, bai_file_id, is_paired) in input_file_ids['bams'].iteritems():
            for reference_subset in grouped_references:
                j = job.addChildJobFn(filter_bam, bam_file_id, bai_file_id, reference_subset, is_paired)
                filtered_bam_file_ids[reference_subset].append(j.rv())
    else:
        filtered_bam_file_ids = {}
    if input_file_ids['annotation'] is not None:
        j = job.addChildJobFn(generate_annotation_hints, input_file_ids['annotation'], genome)
        annotation_hints_file_id = j.rv()
    else:
        annotation_hints_file_id = None
    # returns path to filtered GFF output
    return job.addFollowOnJobFn(merge_bams, filtered_bam_file_ids, annotation_hints_file_id, genome).rv()


def filter_bam(job, bam_file_id, bai_file_id, reference_subset, is_paired):
    """
    Slices out a chromosome from a BAM, re-sorts by name, filters the reads, then re-sorts by position.
    filterBam does weird things when piped to stdout, so I don't do that.
    """
    job.fileStore.logToMaster('Name-sorting and filtering {}'.format(bam_file_id), level=logging.INFO)
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


def merge_bams(job, filtered_bam_file_ids, annotation_hints_file_id, genome):
    """
    Takes a dictionary mapping reference chunks to filtered BAMs. For each reference chunk, these BAMs will be
    first concatenated then sorted, then passed off to hint building.
    """
    job.fileStore.logToMaster('Merging BAMs for {}.'.format(genome), level=logging.INFO)
    # filtered_bam_file_ids may be empty, which signifies that for this genome we have no BAMs
    if len(filtered_bam_file_ids) > 0:
        merged_bam_file_ids = {}
        for ref_group, file_ids in filtered_bam_file_ids.iteritems():
            merged_bam_file_ids[ref_group] = job.addChildJobFn(cat_sort_bams, file_ids, ref_group, genome).rv()
    else:
        merged_bam_file_ids = None
    return job.addFollowOnJobFn(build_hints, merged_bam_file_ids, annotation_hints_file_id, genome).rv()


def cat_sort_bams(job, bam_file_ids, ref_group, genome):
    """
    Takes a list of bam file IDs and combines/sorts them.

    TODO: the 4096 file hack below is hacky. Should only be a problem for very fragmented references.
    """
    job.fileStore.logToMaster('{}: Concatenating and sorting bams for references: {}'.format(genome,
                                                                                             ','.join(ref_group)),
                              level=logging.INFO)
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


def build_hints(job, merged_bam_file_ids, annotation_hints_file_id, genome):
    """
    Takes the merged BAM for a genome and produces both intron and exon hints.
    """
    job.fileStore.logToMaster('Building hints for {}'.format(genome))
    if merged_bam_file_ids is None:
        intron_hints_file_ids = exon_hints_file_ids = None
    else:
        exon_hints_file_ids = {}
        intron_hints_file_ids = {}
        for ref_group, file_ids in merged_bam_file_ids.iteritems():
            intron_hints_file_ids[ref_group] = job.addChildJobFn(build_intron_hints, file_ids).rv()
            exon_hints_file_ids[ref_group] = job.addChildJobFn(build_exon_hints, file_ids).rv()
    return job.addFollowOnJobFn(cat_hints, intron_hints_file_ids, exon_hints_file_ids,
                                annotation_hints_file_id, genome).rv()


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


def generate_annotation_hints(job, annotation_hints_file_id, genome):
    """
    Converts the annotation file into hints. First converts the gff3 directly to genePred so we can make use
    of the transcript library.

    Hints are derived from both CDS exonic intervals and intron intervals
    """
    job.fileStore.logToMaster('Generating annotation hints for {}'.format(genome), level=logging.INFO)
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


def cat_hints(job, intron_hints_file_ids, exon_hints_file_ids, annotation_hints, genome):
    """Returns file ID to combined, sorted hints"""
    if all(x is None for x in (intron_hints_file_ids, exon_hints_file_ids, annotation_hints)):
        return None  # we are just loading the genome
    # load every hint file to this job
    job.fileStore.logToMaster('Concatenating final hints for {}'.format(genome))
    if intron_hints_file_ids is not None:
        intron_hints = [job.fileStore.readGlobalFile(x) for x in intron_hints_file_ids.itervalues()]
    else:
        intron_hints = []
    if exon_hints_file_ids is not None:
        exon_hints = [job.fileStore.readGlobalFile(x) for x in exon_hints_file_ids.itervalues()]
    else:
        exon_hints = []
    if annotation_hints is not None:
        annotation_hints = [job.fileStore.readGlobalFile(annotation_hints)]
    else:
        annotation_hints = []
    # generate a temporary file we will cat everything to.
    cat_hints = tools.fileOps.get_tmp_toil_file()
    with open(cat_hints, 'w') as outf:
        for file in itertools.chain(intron_hints, exon_hints, annotation_hints):
            for line in open(file):
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
    Make sure that this BAM is actually aligned to this fasta
    """
    handle = pysam.Samfile(bam_path, 'rb')
    if set(handle.references) != fasta_sequences:
        base_err = 'Error: BAM {} does not have the same sequences as the FASTA for genome {}'
        err = base_err.format(bam_path, genome)
        raise UserException(err)


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


def group_references(sam_handle, num_bases=10 ** 6, max_seqs=100):
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
