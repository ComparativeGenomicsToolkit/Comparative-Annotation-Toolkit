"""
Comparative Annotation Toolkit.
"""
import string
import random
import datetime
import collections
import itertools
import logging
import os
import shutil
import json
from collections import OrderedDict
from frozendict import frozendict
from configobj import ConfigObj
from subprocess import check_call

import luigi
import luigi.contrib.sqla
from luigi.util import requires
from toil.job import Job
import pandas as pd
from bx.intervals.cluster import ClusterTree

import tools.bio
import tools.fileOps
import tools.intervals
import tools.hal
import tools.misc
import tools.nameConversions
import tools.procOps
import tools.mathOps
import tools.psl
import tools.sqlInterface
import tools.sqlite
import tools.hintsDatabaseInterface
import tools.transcripts
import tools.gff3
from tools.luigiAddons import multiple_requires, IndexTarget
from align_transcripts import align_transcripts
from augustus import augustus
from augustus_cgp import augustus_cgp
from augustus_pb import augustus_pb
from chaining import chaining
from classify import classify
from consensus import generate_consensus, load_alt_names, load_hgm_vectors
from filter_transmap import filter_transmap
from hgm import hgm, parse_hgm_gtf
from transmap_classify import transmap_classify
from plots import generate_plots
from hints_db import hints_db
from parent_gene_assignment import assign_parents
from exceptions import *

logger = logging.getLogger('cat')


###
# Base tasks shared by pipeline tasks
###


class PipelineTask(luigi.Task):
    """
    Base class for all tasks in this pipeline. Has Parameters for all input parameters that will be inherited
    by all downstream tools.

    Provides useful methods for handling parameters being passed between modules.

    Note: significant here is not the same as significant in get_pipeline_args. Significant here is for the luigi
    scheduler to know which parameters define a unique task ID. This would come into play if multiple instances of this
    pipeline are being run on the same scheduler at once.
    """
    hal = luigi.Parameter()
    ref_genome = luigi.Parameter()
    config = luigi.Parameter()
    out_dir = luigi.Parameter(default='./cat_output')
    work_dir = luigi.Parameter(default='./cat_work')
    target_genomes = luigi.TupleParameter(default=None)
    annotate_ancestors = luigi.BoolParameter(default=False)
    binary_mode = luigi.ChoiceParameter(choices=["docker", "local", "singularity"], default='docker',
                                        significant=False)
    # AugustusTM(R) parameters
    augustus = luigi.BoolParameter(default=False)
    augustus_species = luigi.Parameter(default='human', significant=False)
    tm_cfg = luigi.Parameter(default='augustus_cfgs/extrinsic.ETM1.cfg', significant=False)
    tmr_cfg = luigi.Parameter(default='augustus_cfgs/extrinsic.ETM2.cfg', significant=False)
    augustus_utr_off = luigi.BoolParameter(default=False, significant=False)
    # AugustusCGP parameters
    augustus_cgp = luigi.BoolParameter(default=False)
    cgp_param = luigi.Parameter(default=None, significant=False)
    augustus_cgp_cfg_template = luigi.Parameter(default='augustus_cfgs/cgp_extrinsic_template.cfg', significant=False)
    maf_chunksize = luigi.IntParameter(default=2500000, significant=False)
    maf_overlap = luigi.IntParameter(default=500000, significant=False)
    cgp_train_num_exons = luigi.IntParameter(default=5000, significant=False)
    # AugustusPB parameters
    augustus_pb = luigi.BoolParameter(default=False)
    pb_genome_chunksize = luigi.IntParameter(default=5000000, significant=False)
    pb_genome_overlap = luigi.IntParameter(default=500000, significant=False)
    pb_cfg = luigi.Parameter(default='augustus_cfgs/extrinsic.M.RM.PB.E.W.cfg', significant=False)
    # Hgm parameters
    hgm_cpu = luigi.IntParameter(default=4, significant=False)
    # assemblyHub parameters
    assembly_hub = luigi.BoolParameter(default=False)
    hub_email = luigi.Parameter(default='NoEmail', significant=False)
    # Paralogy detection options
    global_near_best = luigi.FloatParameter(default=0.15, significant=False)
    # consensus options
    intron_rnaseq_support = luigi.IntParameter(default=0, significant=False)
    exon_rnaseq_support = luigi.IntParameter(default=0, significant=False)
    intron_annot_support = luigi.IntParameter(default=0, significant=False)
    exon_annot_support = luigi.IntParameter(default=0, significant=False)
    original_intron_support = luigi.IntParameter(default=0, significant=False)
    denovo_num_introns = luigi.IntParameter(default=0, significant=False)
    denovo_splice_support = luigi.IntParameter(default=0, significant=False)
    denovo_exon_support = luigi.IntParameter(default=0, significant=False)
    require_pacbio_support = luigi.BoolParameter(default=False, significant=False)
    in_species_rna_support_only = luigi.BoolParameter(default=False, significant=True)
    rebuild_consensus = luigi.BoolParameter(default=False, significant=True)
    # Toil options
    batchSystem = luigi.Parameter(default='singleMachine', significant=False)
    maxCores = luigi.IntParameter(default=8, significant=False)
    parasolCommand = luigi.Parameter(default=None, significant=False)
    defaultMemory = luigi.Parameter(default='8G', significant=False)
    disableCaching = luigi.BoolParameter(default=False, significant=False)
    workDir = luigi.Parameter(default=None, significant=False)
    defaultDisk = luigi.Parameter(default='8G', significant=False)
    cleanWorkDir = luigi.Parameter(default='onSuccess', significant=False)
    provisioner = luigi.Parameter(default=None, significant=False)
    nodeTypes = luigi.Parameter(default=None, significant=False)
    maxNodes = luigi.Parameter(default=None, significant=False)
    minNode = luigi.Parameter(default=None, significant=False)
    metrices = luigi.Parameter(default=None, significant=False)
    zone = luigi.Parameter(default=None, significant=False)

    def __repr__(self):
        """override the repr to make logging cleaner"""
        if hasattr(self, 'genome'):
            return 'Task: {} for {}'.format(self.__class__.__name__, self.genome)
        elif hasattr(self, 'mode'):
            return 'Task: {} for {}'.format(self.__class__.__name__, self.mode)
        else:
            return 'Task: {}'.format(self.__class__.__name__)

    def get_pipeline_args(self):
        """returns a namespace of all of the arguments to the pipeline. Resolves the target genomes variable"""
        # We use this environment variable as a bit of global state,
        # to avoid threading this through in each of the hundreds of
        # command invocations.
        os.environ['CAT_BINARY_MODE'] = self.binary_mode

        args = tools.misc.PipelineNamespace()
        args.set('binary_mode', self.binary_mode, False)
        args.set('hal', os.path.abspath(self.hal), True)
        args.set('ref_genome', self.ref_genome, True)
        args.set('out_dir', os.path.abspath(self.out_dir), True)
        args.set('work_dir', os.path.abspath(self.work_dir), True)
        args.set('augustus', self.augustus, True)
        args.set('augustus_cgp', self.augustus_cgp, True)
        args.set('augustus_pb', self.augustus_pb, True)
        args.set('augustus_species', self.augustus_species, True)
        args.set('tm_cfg', os.path.abspath(self.tm_cfg), True)
        args.set('tmr_cfg', os.path.abspath(self.tmr_cfg), True)
        args.set('augustus_cgp', self.augustus_cgp, True)
        args.set('maf_chunksize', self.maf_chunksize, True)
        args.set('maf_overlap', self.maf_overlap, True)
        args.set('pb_genome_chunksize', self.pb_genome_chunksize, True)
        args.set('pb_genome_overlap', self.pb_genome_overlap, True)
        args.set('pb_cfg', os.path.abspath(self.pb_cfg), True)

        args.set('augustus_cgp_cfg_template', os.path.abspath(self.augustus_cgp_cfg_template), True)
        args.set('augustus_utr_off', self.augustus_utr_off, True)
        if self.cgp_param is not None:
            args.set('cgp_param', os.path.abspath(self.cgp_param), True)
        else:
            args.set('cgp_param', None, True)
        args.set('cgp_train_num_exons', self.cgp_train_num_exons, True)
        args.set('hgm_cpu', self.hgm_cpu, False)

        # user flags for paralog resolution
        args.set('global_near_best', self.global_near_best, True)
        
        # user specified flags for consensus finding
        args.set('intron_rnaseq_support', self.intron_rnaseq_support, False)
        args.set('exon_rnaseq_support', self.exon_rnaseq_support, False)
        args.set('intron_annot_support', self.intron_annot_support, False)
        args.set('exon_annot_support', self.exon_annot_support, False)
        args.set('original_intron_support', self.original_intron_support, False)
        args.set('denovo_num_introns', self.denovo_num_introns, False)
        args.set('denovo_splice_support', self.denovo_splice_support, False)
        args.set('denovo_exon_support', self.denovo_exon_support, False)
        args.set('require_pacbio_support', self.require_pacbio_support, False)
        args.set('in_species_rna_support_only', self.in_species_rna_support_only, False)
        args.set('rebuild_consensus', self.rebuild_consensus, False)

        # stats location
        args.set('stats_db', os.path.join(args.out_dir, 'databases', 'timing_stats.db'), False)

        # flags for assembly hub building
        args.set('assembly_hub', self.assembly_hub, False)  # assembly hub doesn't need to cause rebuild of gene sets
        args.set('hub_email', self.hub_email, False)

        # flags for figuring out which genomes we are going to annotate
        args.set('annotate_ancestors', self.annotate_ancestors, True)

        # halStats is run below, before any validate() methods are called.
        if not tools.misc.is_exec('halStats'):
            raise ToolMissingException('halStats from the HAL tools package not in global path')

        args.set('hal_genomes', tools.hal.extract_genomes(args.hal, self.annotate_ancestors), True)
        target_genomes = tools.hal.extract_genomes(args.hal, self.annotate_ancestors, self.target_genomes)
        target_genomes = tuple(x for x in target_genomes if x != self.ref_genome)
        args.set('target_genomes', target_genomes, True)

        args.set('cfg', self.parse_cfg(), True)
        args.set('dbs', PipelineTask.get_databases(args), True)
        args.set('annotation', args.cfg['ANNOTATION'][args.ref_genome], True)
        args.set('hints_db', os.path.join(args.work_dir, 'hints_database', 'hints.db'), True)
        args.set('rnaseq_genomes', frozenset(set(args.cfg['INTRONBAM'].keys()) | set(args.cfg['BAM'].keys())), True)
        args.set('intron_only_genomes', frozenset(set(args.cfg['INTRONBAM'].keys()) - set(args.cfg['BAM'].keys())), True)
        args.set('isoseq_genomes', frozenset(args.cfg['ISO_SEQ_BAM'].keys()), True)
        args.set('annotation_genomes', frozenset(args.cfg['ANNOTATION'].keys()), True)
        args.set('modes', self.get_modes(args), True)
        args.set('augustus_tmr', True if 'augTMR' in args.modes else False, True)

        if self.__class__.__name__ in ['RunCat', 'Augustus', 'AugustusCgp', 'AugustusPb']:
            self.validate_cfg(args)

        return args

    def parse_cfg(self):
        """
        Parses the input config file. Config file format:

        [ANNOTATION]
        annotation = /path/to/gff3

        [INTRONBAM]
        genome1 = /path/to/non_polyA_bam1.bam, /path/to/non_polyA_bam2.bam

        [BAM]
        genome1 = /path/to/fofn

        [ISO_SEQ_BAM]
        genome1 = /path/to/bam/or/fofn

        [PROTEIN_FASTA]
        genome1 = /path/to/fasta/or/fofn

        The annotation field must be populated for the reference genome.

        The protein fasta field should be populated for every genome you wish to perform protein alignment to.

        BAM annotations can be put either under INTRONBAM or BAM. Any INTRONBAM will only have intron data loaded,
        and is suitable for lower quality RNA-seq.

        """
        if not os.path.exists(self.config):
            raise MissingFileException('Config file {} not found.'.format(self.config))
        # configspec validates the input config file
        configspec = ['[ANNOTATION]', '__many__ = string',
                      '[INTRONBAM]', '__many__ = list',
                      '[BAM]', '__many__ = list',
                      '[ISO_SEQ_BAM]', '__many__ = list',
                      '[PROTEIN_FASTA]', '__many__ = list']
        parser = ConfigObj(self.config, configspec=configspec)

        # convert the config into a new dict, parsing the fofns
        cfg = collections.defaultdict(dict)
        for dtype in ['ANNOTATION', 'PROTEIN_FASTA']:
            if dtype not in parser:
                cfg[dtype] = {}
            else:
                for genome, annot in parser[dtype].iteritems():
                    annot = os.path.abspath(annot)
                    if not os.path.exists(annot):
                        raise MissingFileException('Missing {} file {}.'.format(dtype.lower(), annot))
                    cfg[dtype][genome] = annot

        # if a given genome only has one BAM, it is a string. Fix this. Extract all paths from fofn files.
        for dtype in ['BAM', 'INTRONBAM', 'ISO_SEQ_BAM']:
            if dtype not in parser:  # the user does not have to specify all field types
                cfg[dtype] = {}
                continue
            for genome in parser[dtype]:
                path = parser[dtype][genome]
                if isinstance(path, str):
                    if not tools.misc.is_bam(path):
                        # this is a fofn
                        cfg[dtype][genome] = [os.path.abspath(x.rstrip()) for x in open(path)]
                    else:
                        # this is a single BAM
                        cfg[dtype][genome] = [os.path.abspath(path)]
                else:
                    cfg[dtype][genome] = []
                    for p in path:
                        if tools.misc.is_bam(p):  # is a bam
                            cfg[dtype][genome].append(os.path.abspath(p))
                        else:  # assume to be a fofn
                            cfg[dtype][genome].extend([os.path.abspath(x.rstrip()) for x in open(p)])

        # return a hashable version
        return frozendict((key, frozendict((ikey, tuple(ival) if isinstance(ival, list) else ival)
                                           for ikey, ival in val.iteritems())) for key, val in cfg.iteritems())

    def validate_cfg(self, args):
        """Validate the input config file."""
        if len(args.cfg['BAM']) + len(args.cfg['INTRONBAM']) + \
                len(args.cfg['ISO_SEQ_BAM']) + len(args.cfg['ANNOTATION']) == 0:
            logger.warning('No extrinsic data or annotations found in config. Will load genomes only.')
        elif len(args.cfg['BAM']) + len(args.cfg['INTRONBAM']) + len(args.cfg['ISO_SEQ_BAM']) == 0:
            logger.warning('No extrinsic data found in config. Will load genomes and annotation only.')

        for dtype in ['BAM', 'INTRONBAM', 'ISO_SEQ_BAM']:
            for genome in args.cfg[dtype]:
                for bam in args.cfg[dtype][genome]:
                    if not os.path.exists(bam):
                        raise MissingFileException('Missing BAM {}.'.format(bam))
                    if not tools.misc.is_bam(bam):
                        raise InvalidInputException('BAM {} is not a valid BAM.'.format(bam))
                    if not os.path.exists(bam + '.bai'):
                        raise MissingFileException('Missing BAM index {}.'.format(bam + '.bai'))

        for dtype in ['ANNOTATION', 'PROTEIN_FASTA']:
            for genome, annot in args.cfg[dtype].iteritems():
                if not os.path.exists(annot):
                    raise MissingFileException('Missing {} file {}.'.format(dtype.lower(), annot))

        if all(g in args.hal_genomes for g in args.target_genomes) is False:
            bad_genomes = set(args.hal_genomes) - set(args.target_genomes)
            err_msg = 'Genomes {} present in configuration and not present in HAL.'.format(','.join(bad_genomes))
            raise UserException(err_msg)

        if args.ref_genome not in args.cfg['ANNOTATION']:
            raise UserException('Reference genome {} did not have a provided annotation.'.format(self.ref_genome))

        # raise if the user if the user is providing dubious inputs
        if args.augustus_cgp and len(args.rnaseq_genomes) == 0:
            raise InvalidInputException('AugustusCGP is being ran without any RNA-seq hints!')
        if args.augustus_pb and len(args.isoseq_genomes) == 0:
            raise InvalidInputException('AugustusPB is being ran without any IsoSeq hints!')

    def get_modes(self, args):
        """returns a tuple of the execution modes being used here"""
        modes = ['transMap']
        if args.augustus_cgp is True:
            modes.append('augCGP')
        if args.augustus is True:
            modes.append('augTM')
            if len(set(args.rnaseq_genomes) & set(args.target_genomes)) > 0:
                modes.append('augTMR')
        if args.augustus_pb is True:
            modes.append('augPB')
        return tuple(modes)

    def get_module_args(self, module, **args):
        """
        convenience wrapper that takes a parent module and propagates any required arguments to generate the full
        argument set.
        """
        pipeline_args = self.get_pipeline_args()
        return module.get_args(pipeline_args, **args)

    @staticmethod
    def get_databases(pipeline_args):
        """wrapper for get_database() that provides all of the databases"""
        dbs = {genome: PipelineTask.get_database(pipeline_args, genome) for genome in pipeline_args.hal_genomes}
        return frozendict(dbs)

    @staticmethod
    def get_database(pipeline_args, genome):
        """database paths must be resolved here to handle multiple programs accessing them"""
        base_out_dir = os.path.join(pipeline_args.out_dir, 'databases')
        return os.path.join(base_out_dir, '{}.db'.format(genome))

    @staticmethod
    def get_plot_dir(pipeline_args, genome):
        """plot base directories must be resolved here to handle multiple programs accessing them"""
        base_out_dir = os.path.join(pipeline_args.out_dir, 'plots')
        return os.path.join(base_out_dir, genome)

    @staticmethod
    def get_metrics_dir(pipeline_args, genome):
        """plot data directories must be resolved here to handle multiple programs accessing them"""
        base_out_dir = os.path.join(pipeline_args.work_dir, 'plot_data')
        return os.path.join(base_out_dir, genome)

    @staticmethod
    def write_metrics(metrics_dict, out_target):
        """write out a metrics dictionary to a path for later loading and plotting"""
        tools.fileOps.ensure_file_dir(out_target.path)
        with out_target.open('w') as outf:
            json.dump(metrics_dict, outf)


@PipelineTask.event_handler(luigi.Event.PROCESSING_TIME)
def processing_time(task, processing_time):
    """
    An event to record processing time of each task. This event records directly to a sqlite database.
    """
    pipeline_args = task.get_pipeline_args()
    stats_db = pipeline_args.stats_db
    finish_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with tools.sqlite.ExclusiveSqlConnection(stats_db) as engine:
        c = engine.cursor()
        c.execute('create table if not exists stats '
                  '(TaskId string unique, FinishTime string, ProcessingTime real)')
        c.execute('insert or replace into stats values (?, ?, ?)', [task.task_id, finish_time, processing_time])
        engine.commit()


class PipelineWrapperTask(PipelineTask, luigi.WrapperTask):
    """add WrapperTask functionality to PipelineTask"""
    pass


class AbstractAtomicFileTask(PipelineTask):
    """
    Abstract Task for single files.
    """
    def run_cmd(self, cmd):
        """
        Run a external command that will produce the output file for this task to stdout. Capture this to the file.
        """
        # luigi localTargets guarantee atomicity if used as a context manager
        with self.output().open('w') as outf:
            tools.procOps.run_proc(cmd, stdout=outf)


class ToilTask(PipelineTask):
    """
    Task for launching toil pipelines from within luigi.
    """
    resources = {'toil': 1}  # all toil pipelines use 1 toil

    def __repr__(self):
        """override the PipelineTask repr to report the batch system being used"""
        base_repr = super(ToilTask, self).__repr__()
        return 'Toil' + base_repr + ' using batchSystem {}'.format(self.batchSystem)

    def prepare_toil_options(self, work_dir):
        """
        Prepares a Namespace object for Toil which has all defaults, overridden as specified
        Will see if the jobStore path exists, and if it does, assume that we need to add --restart
        :param work_dir: Parent directory where toil work will be done. jobStore will be placed inside. Will be used
        to fill in the workDir class variable.
        :return: Namespace
        """
        toil_args = self.get_toil_defaults()
        toil_args.__dict__.update(vars(self))
        toil_args.stats = True
        toil_args.defaultPreemptable = True
        if self.zone is not None:
            #job_store = self.provisioner + ':' + self.zone + ':' + ''.join(random.choice(string.ascii_lowercase) for m in range(7))

            job_dir = os.path.join(work_dir, 'jobStore') # Directory where the AWS directory file is
            if os.path.exists(job_dir):
                for i in os.listdir(job_dir):
                    if os.path.isfile(os.path.join(job_dir,i)) and self.provisioner in i:
		        job_store = i
                        toil_args.restart = True
                        break
            if toil_args.restart is not True:
                job_store = self.provisioner + ':' + self.zone + ':' + ''.join(random.choice(string.ascii_lowercase) for m in range(7))
                try:
                    os.makedirs(job_dir)
                except OSError:
                    pass
                open(os.path.join(job_dir,job_store),'w+').close() # Creates empty file with the title as the AWS jobstore

        else:
            job_store = os.path.join(work_dir, 'jobStore')
            tools.fileOps.ensure_file_dir(job_store)

            # this logic tries to determine if we should try and restart an existing jobStore
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

        if tools.misc.running_in_container():
            # Caching doesn't work in containers, because the
            # container filesystems are transient overlays that don't
            # support hardlinking.
            toil_args.disableCaching = True
        if toil_args.batchSystem == 'parasol' and toil_args.disableCaching is False:
            raise RuntimeError('Running parasol without disabled caching is a very bad idea.')
        if toil_args.batchSystem == 'parasol' and toil_args.workDir is None:
            raise RuntimeError('Running parasol without setting a shared work directory will not work. Please specify '
                               '--workDir.')
        if toil_args.workDir is not None:
            tools.fileOps.ensure_dir(toil_args.workDir)
        #job_store = 'file:' + job_store
        toil_args.jobStore = job_store
        self.job_store = job_store
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


@ToilTask.event_handler(luigi.Event.SUCCESS)
def success(task):
    """
    An event to record the total CPU time of a toil job.
    """
    pipeline_args = task.get_pipeline_args()
    stats_db = pipeline_args.stats_db
    if task.zone is not None:
        cmd = ['toil', 'stats', '--raw', task.job_store]
        try:
            os.remove(os.path.abspath(task.job_store))
        except OSError:
            pass
    else: 
        cmd = ['toil', 'stats', '--raw', os.path.abspath(task.job_store)]
    raw = tools.procOps.call_proc(cmd)
    parsed = raw[raw.index('{'):raw.rfind('}') + 1]
    stats = json.loads(parsed)
    with tools.sqlite.ExclusiveSqlConnection(stats_db) as engine:
        c = engine.cursor()
        c.execute('create table if not exists toil_stats '
                  '(TaskId string unique, TotalTime real, AverageTime real)')
        c.execute('insert or replace into toil_stats values (?, ?, ?)', [task.task_id,
                                                                         stats['jobs']['total_clock'],
                                                                         stats['jobs']['average_clock']])
        engine.commit()


class RebuildableTask(PipelineTask):
    def __init__(self, *args, **kwargs):
        """Allows us to force a task to be re-run. https://github.com/spotify/luigi/issues/595"""
        super(PipelineTask, self).__init__(*args, **kwargs)
        # To force execution, we just remove all outputs before `complete()` is called
        if self.rebuild_consensus is True:
            outputs = luigi.task.flatten(self.output())
            for out in outputs:
                if out.exists():
                    out.remove()


class TrackTask(RebuildableTask):
    """Provides shared values for all of the track tasks"""
    genome = luigi.Parameter()
    track_path = luigi.Parameter()
    trackdb_path = luigi.Parameter()

    def requires(self):
        yield self.clone(Consensus)
        yield self.clone(EvaluateTransMap)
        yield self.clone(EvaluateTranscripts)
        yield self.clone(Hgm)
        yield self.clone(ReferenceFiles)
        yield self.clone(EvaluateTransMap)
        yield self.clone(TransMap)

    def output(self):
        return luigi.LocalTarget(self.track_path), luigi.LocalTarget(self.trackdb_path)


###
# pipeline tasks
###


class RunCat(PipelineWrapperTask):
    """
    Task that executes the entire pipeline.
    """
    def validate(self, pipeline_args):
        """General input validation"""
        if pipeline_args.binary_mode == 'docker':
            # Update docker container
            check_call(['docker', 'pull', 'quay.io/ucsc_cgl/cat:latest'])
        if not os.path.exists(pipeline_args.hal):
            raise InputMissingException('HAL file not found at {}.'.format(pipeline_args.hal))
        for d in [pipeline_args.out_dir, pipeline_args.work_dir]:
            if not os.path.exists(d):
                if not tools.fileOps.dir_is_writeable(os.path.dirname(d)):
                    raise UserException('Cannot create directory {}.'.format(d))
            else:
                if not tools.fileOps.dir_is_writeable(d):
                    raise UserException('Directory {} is not writeable.'.format(d))
        if not os.path.exists(pipeline_args.annotation):
            raise InputMissingException('Annotation file {} not found.'.format(pipeline_args.annotation))
        # TODO: validate augustus species, tm/tmr/cgp/param files.
        if pipeline_args.ref_genome not in pipeline_args.hal_genomes:
            raise InvalidInputException('Reference genome {} not present in HAL.'.format(pipeline_args.ref_genome))
        missing_genomes = {g for g in pipeline_args.target_genomes if g not in pipeline_args.hal_genomes}
        if len(missing_genomes) > 0:
            missing_genomes = ','.join(missing_genomes)
            raise InvalidInputException('Target genomes {} not present in HAL.'.format(missing_genomes))
        if pipeline_args.ref_genome in pipeline_args.target_genomes:
            raise InvalidInputException('A target genome cannot be the reference genome.')

    def requires(self):
        pipeline_args = self.get_pipeline_args()
        self.validate(pipeline_args)
        yield self.clone(PrepareFiles)
        yield self.clone(BuildDb)
        yield self.clone(Chaining)
        yield self.clone(TransMap)
        yield self.clone(EvaluateTransMap)
        if self.augustus is True:
            yield self.clone(Augustus)
        if self.augustus_cgp is True:
            yield self.clone(AugustusCgp)
            yield self.clone(FindDenovoParents, mode='augCGP')
        if self.augustus_pb is True:
            yield self.clone(AugustusPb)
            yield self.clone(FindDenovoParents, mode='augPB')
            yield self.clone(IsoSeqTranscripts)
        yield self.clone(Hgm)
        yield self.clone(AlignTranscripts)
        yield self.clone(EvaluateTranscripts)
        yield self.clone(Consensus)
        yield self.clone(Plots)
        if self.assembly_hub is True:
            yield self.clone(AssemblyHub)
        yield self.clone(ReportStats)


class PrepareFiles(PipelineWrapperTask):
    """
    Wrapper for file preparation tasks GenomeFiles and ReferenceFiles
    """
    def requires(self):
        yield self.clone(GenomeFiles)
        yield self.clone(ReferenceFiles)


class GenomeFiles(PipelineWrapperTask):
    """
    WrapperTask for producing all genome files.

    GenomeFiles -> GenomeFasta -> GenomeTwoBit -> GenomeFlatFasta -> GenomeFastaIndex
                -> GenomeSizes

    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'genome_files')
        args = tools.misc.HashableNamespace()
        args.genome = genome
        args.fasta = os.path.join(base_dir, genome + '.fa')
        args.two_bit = os.path.join(base_dir, genome + '.2bit')
        args.sizes = os.path.join(base_dir, genome + '.chrom.sizes')
        args.flat_fasta = os.path.join(base_dir, genome + '.fa.flat')
        return args

    def validate(self):
        for haltool in ['hal2fasta', 'halStats']:
            if not tools.misc.is_exec(haltool):
                    raise ToolMissingException('{} from the HAL tools package not in global path'.format(haltool))
        if not tools.misc.is_exec('faToTwoBit'):
            raise ToolMissingException('faToTwoBit tool from the Kent tools package not in global path.')
        if not tools.misc.is_exec('pyfasta'):
            raise ToolMissingException('pyfasta wrapper not found in global path.')

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for genome in list(pipeline_args.target_genomes) + [pipeline_args.ref_genome]:
            args = self.get_args(pipeline_args, genome)
            yield self.clone(GenomeFasta, **vars(args))
            yield self.clone(GenomeTwoBit, **vars(args))
            yield self.clone(GenomeSizes, **vars(args))
            yield self.clone(GenomeFlatFasta, **vars(args))


class GenomeFasta(AbstractAtomicFileTask):
    """
    Produce a fasta file from a hal file. Requires hal2fasta.
    """
    genome = luigi.Parameter()
    fasta = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.fasta)

    def run(self):
        logger.info('Extracting fasta for {}.'.format(self.genome))
        cmd = ['hal2fasta', os.path.abspath(self.hal), self.genome]
        self.run_cmd(cmd)


@requires(GenomeFasta)
class GenomeTwoBit(AbstractAtomicFileTask):
    """
    Produce a 2bit file from a fasta file. Requires kent tool faToTwoBit.
    Needs to be done BEFORE we flatten.
    """
    two_bit = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.two_bit)

    def run(self):
        logger.info('Converting fasta for {} to 2bit.'.format(self.genome))
        cmd = ['faToTwoBit', self.fasta, '/dev/stdout']
        self.run_cmd(cmd)


class GenomeSizes(AbstractAtomicFileTask):
    """
    Produces a genome chromosome sizes file. Requires halStats.
    """
    genome = luigi.Parameter()
    sizes = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.sizes)

    def run(self):
        logger.info('Extracting chromosome sizes for {}.'.format(self.genome))
        cmd = ['halStats', '--chromSizes', self.genome, os.path.abspath(self.hal)]
        self.run_cmd(cmd)


@requires(GenomeTwoBit)
class GenomeFlatFasta(AbstractAtomicFileTask):
    """
    Flattens a genome fasta in-place using pyfasta. Requires the pyfasta package.
    """
    flat_fasta = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.flat_fasta)

    def run(self):
        logger.info('Flattening fasta for {}.'.format(self.genome))
        cmd = ['pyfasta', 'flatten', self.fasta]
        tools.procOps.run_proc(cmd)


class ReferenceFiles(PipelineWrapperTask):
    """
    WrapperTask for producing annotation files.

    ReferenceFiles -> Gff3ToGenePred -> TranscriptBed -> TranscriptFasta -> FlatTranscriptFasta
                            V
                         FakePsl, TranscriptGtf
    """
    @staticmethod
    def get_args(pipeline_args):
        base_dir = os.path.join(pipeline_args.work_dir, 'reference')
        annotation = os.path.splitext(os.path.basename(pipeline_args.annotation))[0]
        args = tools.misc.HashableNamespace()
        args.annotation_gp = os.path.join(base_dir, annotation + '.gp')
        args.annotation_attrs = os.path.join(base_dir, annotation + '.gp_attrs')
        args.annotation_gtf = os.path.join(base_dir, annotation + '.gtf')
        args.transcript_fasta = os.path.join(base_dir, annotation + '.fa')
        args.transcript_flat_fasta = os.path.join(base_dir, annotation + '.fa.flat')
        args.transcript_bed = os.path.join(base_dir, annotation + '.bed')
        args.duplicates = os.path.join(base_dir, annotation + '.duplicates.txt')
        args.ref_psl = os.path.join(base_dir, annotation + '.psl')
        args.__dict__.update(**vars(GenomeFiles.get_args(pipeline_args, pipeline_args.ref_genome)))
        return args

    def validate(self):
        for tool in ['gff3ToGenePred', 'genePredToBed', 'genePredToFakePsl']:
            if not tools.misc.is_exec(tool):
                    raise ToolMissingException('{} from the Kent tools package not in global path'.format(tool))

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        args = self.get_args(pipeline_args)
        yield self.clone(Gff3ToGenePred, **vars(args))
        yield self.clone(Gff3ToAttrs, **vars(args))
        yield self.clone(TranscriptBed, **vars(args))
        yield self.clone(TranscriptFasta, **vars(args))
        yield self.clone(TranscriptGtf, **vars(args))
        yield self.clone(FlatTranscriptFasta, **vars(args))
        yield self.clone(FakePsl, **vars(args))


class Gff3ToGenePred(AbstractAtomicFileTask):
    """
    Generates a genePred from a gff3 file.
    """
    annotation_gp = luigi.Parameter()
    annotation_attrs = luigi.Parameter()
    duplicates = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.annotation_gp)

    def validate(self):
        c = collections.Counter()
        for l in open(self.output().path):
            l = l.split()
            c[l[0]] += 1
        duplicates = {x for x, y in c.iteritems() if y > 1}
        if len(duplicates) > 0:
            with open(self.duplicates, 'w') as outf:
                for l in duplicates:
                    outf.write(l + '\n')
            raise InvalidInputException('Found {:,} duplicate transcript IDs after parsing input GFF3. '
                                        'Please check your input. One possible cause is the lack of a transcript-level '
                                        'identifier on a gene record. Duplicate IDs have been written to: '
                                        '{}'.format(len(duplicates), self.duplicates))

    def run(self):
        pipeline_args = self.get_pipeline_args()
        logger.info('Converting annotation gff3 to genePred.')
        cmd = tools.gff3.convert_gff3_cmd(self.annotation_attrs, pipeline_args.annotation)
        self.run_cmd(cmd)
        self.validate()


@requires(Gff3ToGenePred)
class Gff3ToAttrs(PipelineTask):
    """
    Converts the attrs file from -attrsOut in gff3ToGenePred into a SQLite table.
    """
    table = tools.sqlInterface.Annotation.__tablename__

    def output(self):
        pipeline_args = self.get_pipeline_args()
        database = pipeline_args.dbs[pipeline_args.ref_genome]
        tools.fileOps.ensure_file_dir(database)
        conn_str = 'sqlite:///{}'.format(database)
        digest = tools.fileOps.hashfile(pipeline_args.annotation)
        attrs_table = luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                          target_table=self.table,
                                                          update_id='_'.join([self.table, digest]))
        return attrs_table

    def requires(self):
        pipeline_args = self.get_pipeline_args()
        return self.clone(Gff3ToGenePred, annotation_gp=ReferenceFiles.get_args(pipeline_args).annotation_gp)

    def run(self):
        logger.info('Extracting gff3 attributes to sqlite database.')
        pipeline_args = self.get_pipeline_args()
        df = tools.gff3.parse_gff3(self.annotation_attrs, self.annotation_gp)
        if 'protein_coding' not in set(df.GeneBiotype) or 'protein_coding' not in set(df.TranscriptBiotype):
            if pipeline_args.augustus:
                raise InvalidInputException('No protein_coding annotations found. This will cause problems for '
                                            'AugustusTMR. Please check your GFF3 input.')
            else:
                logger.critical('No protein_coding annotations found!')
        # validate number parsed
        tot_genes = len(open(self.annotation_gp).readlines())
        if tot_genes != len(df):
            raise InvalidInputException('The number of genes parsed from the attrs file is not the same number as '
                                        'in the genePred. This is a parser failure. Contact Ian and make him fix it.')
        database = pipeline_args.dbs[pipeline_args.ref_genome]
        with tools.sqlite.ExclusiveSqlConnection(database) as engine:
            df.to_sql(self.table, engine, if_exists='replace')
        self.output().touch()


@requires(Gff3ToGenePred)
class TranscriptBed(AbstractAtomicFileTask):
    """
    Produces a BED record from the input genePred annotation. Makes use of Kent tool genePredToBed
    """
    transcript_bed = luigi.Parameter()
    annotation_gp = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.transcript_bed)

    def run(self):
        logger.info('Converting annotation genePred to BED.')
        cmd = ['genePredToBed', self.annotation_gp, '/dev/stdout']
        self.run_cmd(cmd)


@multiple_requires(GenomeFlatFasta, TranscriptBed)
class TranscriptFasta(AbstractAtomicFileTask):
    """
    Produces a fasta for each transcript.
    """
    transcript_fasta = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.transcript_fasta)

    def run(self):
        logger.info('Extracting reference annotation fasta.')
        seq_dict = tools.bio.get_sequence_dict(self.fasta, upper=False)
        seqs = {tx.name: tx.get_mrna(seq_dict) for tx in tools.transcripts.transcript_iterator(self.transcript_bed)}
        with self.output().open('w') as outf:
            for name, seq in seqs.iteritems():
                tools.bio.write_fasta(outf, name, seq)


@requires(Gff3ToGenePred)
class TranscriptGtf(AbstractAtomicFileTask):
    """
    Produces a GTF out of the genePred for the reference
    """
    annotation_gtf = luigi.Parameter()
    annotation_gp = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.annotation_gtf)

    def run(self):
        logger.info('Extracting reference annotation GTF.')
        tools.misc.convert_gp_gtf(self.output(), luigi.LocalTarget(self.annotation_gp))


@requires(TranscriptFasta)
class FlatTranscriptFasta(AbstractAtomicFileTask):
    """
    Flattens the transcript fasta for pyfasta.
    """
    transcript_fasta = luigi.Parameter()
    transcript_flat_fasta = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.transcript_flat_fasta)

    def run(self):
        logger.info('Flattening reference annotation fasta.')
        cmd = ['pyfasta', 'flatten', self.transcript_fasta]
        tools.procOps.run_proc(cmd)


@multiple_requires(Gff3ToGenePred, GenomeSizes)
class FakePsl(AbstractAtomicFileTask):
    """
    Produces a fake PSL mapping transcripts to the genome, using the Kent tool genePredToFakePsl
    """
    ref_psl = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.ref_psl)

    def run(self):
        logger.info('Generating annotation fake PSL.')
        cmd = ['genePredToFakePsl', '-chromSize={}'.format(self.sizes), 'noDB',
               self.annotation_gp, '/dev/stdout', '/dev/null']
        self.run_cmd(cmd)


class BuildDb(PipelineTask):
    """
    Constructs the hints database from a series of reduced hints GFFs.

    TODO: output() should be way smarter than this. Currently, it only checks if the indices have been created.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'hints_database')
        args = tools.misc.HashableNamespace()
        args.genome = genome
        args.fasta = GenomeFiles.get_args(pipeline_args, genome).fasta
        args.hal = pipeline_args.hal
        args.cfg = pipeline_args.cfg
        args.annotation = pipeline_args.cfg['ANNOTATION'].get(genome, None)
        args.protein_fasta = pipeline_args.cfg['PROTEIN_FASTA'].get(genome, None)
        args.hints_path = os.path.join(base_dir, genome + '.extrinsic_hints.gff')
        return args

    def validate(self):
        tools.misc.samtools_version()  # validate samtools version
        for tool in ['load2sqlitedb', 'samtools', 'filterBam', 'bam2hints', 'bam2wig', 'wig2hints.pl', 'bam2hints',
                     'bamToPsl', 'blat2hints.pl', 'gff3ToGenePred', 'join_mult_hints.pl', 'sambamba']:
            if not tools.misc.is_exec(tool):
                raise ToolMissingException('Auxiliary program {} not found on path.'.format(tool))

    def requires(self):
        pipeline_args = self.get_pipeline_args()
        for genome in list(pipeline_args.target_genomes) + [pipeline_args.ref_genome]:
            hints_args = BuildDb.get_args(pipeline_args, genome)
            yield self.clone(GenerateHints, hints_args=hints_args, genome=genome)

    def output(self):
        pipeline_args = self.get_pipeline_args()
        tools.fileOps.ensure_file_dir(pipeline_args.hints_db)
        return IndexTarget(pipeline_args.hints_db)

    def run(self):
        pipeline_args = self.get_pipeline_args()
        self.validate()
        for genome in list(pipeline_args.target_genomes) + [pipeline_args.ref_genome]:
            args = BuildDb.get_args(pipeline_args, genome)
            logger.info('Loading sequence for {} into database.'.format(genome))
            base_cmd = ['load2sqlitedb', '--noIdx', '--clean', '--species={}'.format(genome),
                        '--dbaccess={}'.format(pipeline_args.hints_db)]
            tools.procOps.run_proc(base_cmd + [args.fasta], stdout='/dev/null', stderr='/dev/null')
            if os.path.getsize(args.hints_path) != 0:
                logger.info('Loading hints for {} into database.'.format(genome))
                tools.procOps.run_proc(base_cmd + [args.hints_path], stderr='/dev/null')
        logger.info('Indexing database.')
        cmd = ['load2sqlitedb', '--makeIdx', '--clean', '--dbaccess={}'.format(pipeline_args.hints_db)]
        tools.procOps.run_proc(cmd, stdout='/dev/null', stderr='/dev/null')
        logger.info('Hints database completed.')


class GenerateHints(ToilTask):
    """
    Generate hints for each genome as a separate Toil pipeline.
    """
    hints_args = luigi.Parameter()
    genome = luigi.Parameter()
    stats = luigi.BoolParameter()

    def output(self):
        return luigi.LocalTarget(self.hints_args.hints_path)

    def requires(self):
        return self.clone(PrepareFiles), self.clone(ReferenceFiles)

    def validate(self):
        for tool in ['samtools', 'sambamba']:
            if not tools.misc.is_exec(tool):
                raise ToolMissingException('{} is not in global path.'.format(tool))
        for tool in ['gff3ToGenePred', 'bamToPsl']:
            if not tools.misc.is_exec(tool):
                raise ToolMissingException('{} from the Kent tools package not in global path.'.format(tool))
        for tool in ['join_mult_hints.pl', 'blat2hints.pl', 'wig2hints.pl', 'bam2wig', 'bam2hints', 'filterBam']:
            if not tools.misc.is_exec(tool):
                raise ToolMissingException('{} from the augustus tool package not in global path.'.format(tool))

    def run(self):
        self.validate()
        logger.info('Beginning GenerateHints Toil pipeline for {}.'.format(self.genome))
        work_dir = os.path.abspath(os.path.join(self.work_dir, 'toil', 'hints_db', self.genome))
        toil_options = self.prepare_toil_options(work_dir)
        hints_db(self.hints_args, toil_options)
        logger.info('Finished GenerateHints Toil pipeline for {}.'.format(self.genome))


class Chaining(ToilTask):
    """
    Task that launches the Chaining toil pipeline. This pipeline operates on all genomes at once to reduce the
    repeated downloading of the HAL file.
    """
    @staticmethod
    def get_args(pipeline_args):
        base_dir = os.path.join(pipeline_args.work_dir, 'chaining')
        ref_files = GenomeFiles.get_args(pipeline_args, pipeline_args.ref_genome)
        tgt_files = {genome: GenomeFiles.get_args(pipeline_args, genome) for genome in pipeline_args.target_genomes}
        tgt_two_bits = {genome: tgt_files[genome].two_bit for genome in pipeline_args.target_genomes}
        chain_files = {genome: os.path.join(base_dir, '{}-{}.chain'.format(pipeline_args.ref_genome, genome))
                       for genome in pipeline_args.target_genomes}
        args = tools.misc.HashableNamespace()
        args.hal = pipeline_args.hal
        args.ref_genome = pipeline_args.ref_genome
        args.query_two_bit = ref_files.two_bit
        args.query_sizes = ref_files.sizes
        args.target_two_bits = tgt_two_bits
        args.chain_files = chain_files
        return args

    def output(self):
        pipeline_args = self.get_pipeline_args()
        chain_args = self.get_args(pipeline_args)
        for path in chain_args.chain_files.itervalues():
            yield luigi.LocalTarget(path)

    def validate(self):
        if not tools.misc.is_exec('halLiftover'):
            raise ToolMissingException('halLiftover from the halTools package not in global path.')
        for tool in ['pslPosTarget', 'axtChain', 'chainMergeSort']:
            if not tools.misc.is_exec(tool):
                    raise ToolMissingException('{} from the Kent tools package not in global path.'.format(tool))

    def requires(self):
        yield self.clone(PrepareFiles)

    def run(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        logger.info('Launching Pairwise Chaining toil pipeline.')
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'chaining')
        toil_options = self.prepare_toil_options(toil_work_dir)
        chain_args = self.get_args(pipeline_args)
        chaining(chain_args, toil_options)
        logger.info('Pairwise Chaining toil pipeline is complete.')


class TransMap(PipelineWrapperTask):
    """
    Runs transMap.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'transMap')
        ref_files = ReferenceFiles.get_args(pipeline_args)
        args = tools.misc.HashableNamespace()
        args.two_bit = GenomeFiles.get_args(pipeline_args, genome).two_bit
        args.chain_file = Chaining.get_args(pipeline_args).chain_files[genome]
        args.transcript_fasta = ref_files.transcript_fasta
        args.ref_psl = ref_files.ref_psl
        args.annotation_gp = ref_files.annotation_gp
        args.tm_psl = os.path.join(base_dir, genome + '.psl')
        args.tm_gp = os.path.join(base_dir, genome + '.gp')
        args.tm_gtf = os.path.join(base_dir, genome + '.gtf')
        args.filtered_tm_psl = os.path.join(base_dir, genome + '.filtered.psl')
        args.filtered_tm_gp = os.path.join(base_dir, genome + '.filtered.gp')
        args.metrics_json = os.path.join(PipelineTask.get_metrics_dir(pipeline_args, genome), 'filter_tm_metrics.json')
        args.ref_db_path = pipeline_args.dbs[pipeline_args.ref_genome]
        args.db_path = pipeline_args.dbs[genome]
        args.global_near_best = pipeline_args.global_near_best
        return args

    def validate(self):
        for tool in ['pslMap', 'pslRecalcMatch', 'pslMapPostChain', 'pslCDnaFilter', 'clusterGenes']:
            if not tools.misc.is_exec(tool):
                    raise ToolMissingException('{} from the Kent tools package not in global path.'.format(tool))

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            yield self.clone(TransMapPsl, genome=target_genome)
            yield self.clone(FilterTransMap, genome=target_genome)
            yield self.clone(TransMapGtf, genome=target_genome)


class TransMapPsl(PipelineTask):
    """
    Runs transMap. Requires Kent tools pslMap, pslMapPostChain, pslRecalcMatch
    """
    genome = luigi.Parameter()

    def output(self):
        tm_args = self.get_module_args(TransMap, genome=self.genome)
        return luigi.LocalTarget(tm_args.tm_psl), luigi.LocalTarget(tm_args.tm_gp)

    def requires(self):
        return self.clone(PrepareFiles), self.clone(Chaining), self.clone(ReferenceFiles)

    def run(self):
        tm_args = self.get_module_args(TransMap, genome=self.genome)
        logger.info('Running transMap for {}.'.format(self.genome))
        cmd = [['pslMap', '-chainMapFile', tm_args.ref_psl, tm_args.chain_file, '/dev/stdout'],
               ['pslMapPostChain', '/dev/stdin', '/dev/stdout'],
               ['sort', '-k14,14', '-k16,16n'],
               ['pslRecalcMatch', '/dev/stdin', tm_args.two_bit, tm_args.transcript_fasta, 'stdout'],
               ['sort', '-k10,10']]  # re-sort back to query name for filtering
        tmp_file = luigi.LocalTarget(is_tmp=True)
        with tmp_file.open('w') as tmp_fh:
            tools.procOps.run_proc(cmd, stdout=tmp_fh, stderr='/dev/null')
        tm_psl_tgt, tm_gp_tgt = self.output()
        tools.fileOps.ensure_file_dir(tm_psl_tgt.path)
        with tm_psl_tgt.open('w') as outf:
            for psl_rec in tools.psl.psl_iterator(tmp_file.path, make_unique=True):
                tools.fileOps.print_row(outf, psl_rec.psl_string())
        with tm_gp_tgt.open('w') as outf:
            cmd = ['transMapPslToGenePred', '-nonCodingGapFillMax=80', '-codingGapFillMax=50',
                   tm_args.annotation_gp, tm_psl_tgt.path, '/dev/stdout']
            tools.procOps.run_proc(cmd, stdout=outf)


@requires(TransMapPsl)
class FilterTransMap(PipelineTask):
    """
    Filters transMap output using the localNearBest algorithm.
    """
    eval_table = tools.sqlInterface.TmFilterEval.__tablename__

    def output(self):
        pipeline_args = self.get_pipeline_args()
        tm_args = self.get_module_args(TransMap, genome=self.genome)
        tools.fileOps.ensure_file_dir(tm_args.db_path)
        conn_str = 'sqlite:///{}'.format(tm_args.db_path)
        tm_args = self.get_module_args(TransMap, genome=self.genome)
        return (luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                    target_table=self.eval_table,
                                                    update_id='_'.join([self.eval_table, str(hash(pipeline_args))])),
                luigi.LocalTarget(tm_args.filtered_tm_psl),
                luigi.LocalTarget(tm_args.metrics_json),
                luigi.LocalTarget(tm_args.filtered_tm_gp))

    def run(self):
        tm_args = self.get_module_args(TransMap, genome=self.genome)
        logger.info('Filtering transMap PSL for {}.'.format(self.genome))
        table_target, psl_target, json_target, gp_target = self.output()
        resolved_df = filter_transmap(tm_args.tm_psl, tm_args.ref_psl, tm_args.tm_gp,
                                      tm_args.ref_db_path, psl_target, tm_args.global_near_best, json_target)
        with tools.sqlite.ExclusiveSqlConnection(tm_args.db_path) as engine:
            resolved_df.to_sql(self.eval_table, engine, if_exists='replace')
            table_target.touch()
        with gp_target.open('w') as outf:
            cmd = ['transMapPslToGenePred', '-nonCodingGapFillMax=80', '-codingGapFillMax=50',
                   tm_args.annotation_gp, psl_target.path, '/dev/stdout']
            tools.procOps.run_proc(cmd, stdout=outf)


@requires(FilterTransMap)
class TransMapGtf(PipelineTask):
    """
    Converts the unfiltered transMap PSL to GTF
    """
    def output(self):
        tm_args = self.get_module_args(TransMap, genome=self.genome)
        return luigi.LocalTarget(tm_args.tm_gtf)

    def run(self):
        tm_args = self.get_module_args(TransMap, genome=self.genome)
        logger.info('Creating unfiltered transMap GTF for {}.'.format(self.genome))
        tmp_gp = luigi.LocalTarget(is_tmp=True)
        cmd = ['transMapPslToGenePred', '-nonCodingGapFillMax=80', '-codingGapFillMax=50',
               tm_args.annotation_gp, tm_args.tm_psl, tmp_gp.path]
        tools.procOps.run_proc(cmd)
        tools.misc.convert_gp_gtf(self.output(), tmp_gp)


class EvaluateTransMap(PipelineWrapperTask):
    """
    Evaluates transMap alignments.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        tm_args = TransMap.get_args(pipeline_args, genome)
        args = tools.misc.HashableNamespace()
        args.db_path = pipeline_args.dbs[genome]
        args.filtered_tm_psl = tm_args.filtered_tm_psl
        args.ref_psl = ReferenceFiles.get_args(pipeline_args).ref_psl
        args.filtered_tm_gp = tm_args.filtered_tm_gp
        args.annotation_gp = tm_args.annotation_gp
        args.annotation_gp = ReferenceFiles.get_args(pipeline_args).annotation_gp
        args.genome = genome
        args.fasta = GenomeFiles.get_args(pipeline_args, genome).fasta
        args.ref_genome = pipeline_args.ref_genome
        return args

    def validate(self):
        pass

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            tm_eval_args = EvaluateTransMap.get_args(pipeline_args, target_genome)
            yield self.clone(EvaluateTransMapDriverTask, tm_eval_args=tm_eval_args, genome=target_genome)


class EvaluateTransMapDriverTask(PipelineTask):
    """
    Task for per-genome launching of a toil pipeline for aligning transcripts to their parent.
    """
    genome = luigi.Parameter()
    tm_eval_args = luigi.Parameter()
    table = tools.sqlInterface.TmEval.__tablename__

    def write_to_sql(self, df):
        """Load the results into the SQLite database"""
        with tools.sqlite.ExclusiveSqlConnection(self.tm_eval_args.db_path) as engine:
            df.to_sql(self.table, engine, if_exists='replace')
            self.output().touch()
            logger.info('Loaded table: {}.{}'.format(self.genome, self.table))

    def output(self):
        pipeline_args = self.get_pipeline_args()
        tools.fileOps.ensure_file_dir(self.tm_eval_args.db_path)
        conn_str = 'sqlite:///{}'.format(self.tm_eval_args.db_path)
        return luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                   target_table=self.table,
                                                   update_id='_'.join([self.table, str(hash(pipeline_args))]))

    def requires(self):
        return self.clone(TransMap), self.clone(ReferenceFiles)

    def run(self):
        logger.info('Evaluating transMap results for {}.'.format(self.genome))
        results = transmap_classify(self.tm_eval_args)
        self.write_to_sql(results)


class Augustus(PipelineWrapperTask):
    """
    Runs AugustusTM(R) on the coding output from transMap.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'augustus')
        args = tools.misc.HashableNamespace()
        args.ref_genome = pipeline_args.ref_genome
        args.genome = genome
        args.genome_fasta = GenomeFiles.get_args(pipeline_args, genome).fasta
        args.annotation_gp = ReferenceFiles.get_args(pipeline_args).annotation_gp
        args.filtered_tm_gp = TransMap.get_args(pipeline_args, genome).filtered_tm_gp
        tm_args = TransMap.get_args(pipeline_args, genome)
        args.ref_psl = tm_args.ref_psl
        args.filtered_tm_psl = tm_args.filtered_tm_psl
        args.augustus_tm_gp = os.path.join(base_dir, genome + '.augTM.gp')
        args.augustus_tm_gtf = os.path.join(base_dir, genome + '.augTM.gtf')
        args.tm_cfg = pipeline_args.tm_cfg
        args.tmr_cfg = pipeline_args.tmr_cfg
        args.augustus_species = pipeline_args.augustus_species
        # invert the UTR flag
        args.utr = not pipeline_args.augustus_utr_off
        args.augustus_hints_db = pipeline_args.hints_db
        args.augustus_tmr = genome in pipeline_args.rnaseq_genomes
        if args.augustus_tmr:
            args.augustus_tmr_gp = os.path.join(base_dir, genome + '.augTMR.gp')
            args.augustus_tmr_gtf = os.path.join(base_dir, genome + '.augTMR.gtf')
        return args

    def validate(self):
        for tool in ['augustus', 'transMap2hints.pl']:
            if not tools.misc.is_exec(tool):
                raise ToolMissingException('Auxiliary program {} from the Augustus package not in path.'.format(tool))

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            yield self.clone(AugustusDriverTask, genome=target_genome)


class AugustusDriverTask(ToilTask):
    """
    Task for per-genome launching of a toil pipeline for running Augustus.
    """
    genome = luigi.Parameter()

    def output(self):
        pipeline_args = self.get_pipeline_args()
        augustus_args = Augustus.get_args(pipeline_args, self.genome)
        yield luigi.LocalTarget(augustus_args.augustus_tm_gp)
        yield luigi.LocalTarget(augustus_args.augustus_tm_gtf)
        if augustus_args.augustus_tmr:
            yield luigi.LocalTarget(augustus_args.augustus_tmr_gp)
            yield luigi.LocalTarget(augustus_args.augustus_tmr_gtf)

    def requires(self):
        return self.clone(TransMap), self.clone(BuildDb)

    def extract_coding_genes(self, augustus_args):
        """extracts only coding genes from the input genePred, returning a path to a tmp file"""
        coding_gp = tools.fileOps.get_tmp_file()
        with open(coding_gp, 'w') as outf:
            for tx in tools.transcripts.gene_pred_iterator(augustus_args.filtered_tm_gp):
                if tx.cds_size > 0:
                    tools.fileOps.print_row(outf, tx.get_gene_pred())
        if os.path.getsize(coding_gp) == 0:
            raise InvalidInputException('Unable to extract coding transcripts from the filtered transMap genePred.')
        return coding_gp

    def run(self):
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'augustus', self.genome)
        logger.info('Launching AugustusTMR toil pipeline on {}.'.format(self.genome))
        toil_options = self.prepare_toil_options(toil_work_dir)
        augustus_args = self.get_module_args(Augustus, genome=self.genome)
        coding_gp = self.extract_coding_genes(augustus_args)
        augustus(augustus_args, coding_gp, toil_options)
        logger.info('Augustus toil pipeline for {} completed.'.format(self.genome))
        os.remove(coding_gp)
        for out_gp, out_gtf in tools.misc.pairwise(self.output()):
            tools.misc.convert_gtf_gp(out_gp, out_gtf)


class AugustusCgp(ToilTask):
    """
    Task for launching the AugustusCGP toil pipeline
    """
    @staticmethod
    def get_args(pipeline_args):
        genomes = list(pipeline_args.target_genomes) + [pipeline_args.ref_genome]
        fasta_files = {genome: GenomeFiles.get_args(pipeline_args, genome).fasta for genome in genomes}
        base_dir = os.path.join(pipeline_args.work_dir, 'augustus_cgp')
        # output
        output_gp_files = {genome: os.path.join(base_dir, genome + '.augCGP.gp') for genome in genomes}
        output_gtf_files = {genome: os.path.join(base_dir, genome + '.augCGP.gtf') for genome in genomes}
        raw_output_gtf_files = {genome: os.path.join(base_dir, genome + '.raw.augCGP.gtf') for genome in genomes}
        args = tools.misc.HashableNamespace()
        args.genomes = genomes
        args.annotate_ancestors = pipeline_args.annotate_ancestors
        args.fasta_files = fasta_files
        args.hal = pipeline_args.hal
        args.ref_genome = pipeline_args.ref_genome
        args.augustus_cgp_gp = output_gp_files
        args.augustus_cgp_gtf = output_gtf_files
        args.augustus_cgp_raw_gtf = raw_output_gtf_files
        args.stdout_file = os.path.join(base_dir, 'CGP_stdout.txt')
        args.species = pipeline_args.augustus_species
        args.chunksize = pipeline_args.maf_chunksize
        args.overlap = pipeline_args.maf_overlap
        args.cgp_param = pipeline_args.cgp_param
        if args.cgp_param is None:
            args.param_out_path = os.path.join(base_dir, 'trained_parameters.cfg')
        args.num_exons = pipeline_args.cgp_train_num_exons
        args.hints_db = pipeline_args.hints_db
        args.query_sizes = GenomeFiles.get_args(pipeline_args, pipeline_args.ref_genome).sizes
        args.gtf = ReferenceFiles.get_args(pipeline_args).annotation_gtf
        return args

    def output(self):
        pipeline_args = self.get_pipeline_args()
        cgp_args = self.get_args(pipeline_args)
        for path_dict in [cgp_args.augustus_cgp_gp, cgp_args.augustus_cgp_gtf, cgp_args.augustus_cgp_raw_gtf]:
            for path in path_dict.itervalues():
                yield luigi.LocalTarget(path)

    def validate(self):
        for tool in ['joingenes', 'augustus', 'hal2maf']:
            if not tools.misc.is_exec(tool):
                raise ToolMissingException('tool {} not in global path.'.format(tool))

    def requires(self):
        yield self.clone(TransMap), self.clone(ReferenceFiles), self.clone(BuildDb)

    def prepare_cgp_cfg(self, pipeline_args):
        """use the config template to create a config file"""
        # bam genomes have IsoSeq and/or at least one BAM
        bam_genomes = (pipeline_args.rnaseq_genomes | pipeline_args.isoseq_genomes) - \
                      (pipeline_args.annotation_genomes | pipeline_args.intron_only_genomes)
        # intron only genomes have only intron hints
        intron_only_genomes = pipeline_args.intron_only_genomes - (bam_genomes | pipeline_args.annotation_genomes)
        if not tools.mathOps.all_disjoint([bam_genomes, intron_only_genomes, pipeline_args.annotation_genomes]):
            raise UserException('Error in CGP configuration. Not all genome groups are disjoint.')
        # if --target-genomes is set, remove these genomes from the groups
        target_genomes = set(pipeline_args.target_genomes)
        target_genomes.add(pipeline_args.ref_genome)
        annotation_genomes = pipeline_args.annotation_genomes & target_genomes
        bam_genomes = bam_genomes & target_genomes
        intron_only_genomes = intron_only_genomes & target_genomes
        annotation_genomes = 'none' if len(pipeline_args.annotation_genomes) == 0 else ' '.join(annotation_genomes)
        bam_genomes = 'none' if len(bam_genomes) == 0 else ' '.join(bam_genomes)
        intron_only_genomes = 'none' if len(intron_only_genomes) == 0 else ' '.join(intron_only_genomes)
        template = open(pipeline_args.augustus_cgp_cfg_template).read()
        cfg = template.format(annotation_genomes=annotation_genomes, target_genomes=bam_genomes,
                              intron_target_genomes=intron_only_genomes)
        out_path = tools.fileOps.get_tmp_file()
        with open(out_path, 'w') as outf:
            outf.write(cfg)
        return out_path

    def run(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        logger.info('Launching AugustusCGP toil pipeline.')
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'augustus_cgp')
        toil_options = self.prepare_toil_options(toil_work_dir)
        cgp_args = self.get_args(pipeline_args)
        cgp_args.cgp_cfg = self.prepare_cgp_cfg(pipeline_args)
        augustus_cgp(cgp_args, toil_options)
        logger.info('Finished AugustusCGP toil pipeline.')


class AugustusPb(PipelineWrapperTask):
    """
    Runs AugustusPB. This mode is done on a per-genome basis, but ignores transMap information and and relies only on
    a combination of IsoSeq and RNA-seq
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'augustus_pb')
        args = tools.misc.HashableNamespace()
        args.genome = genome
        genome_files = GenomeFiles.get_args(pipeline_args, genome)
        args.genome_fasta = genome_files.fasta
        args.chrom_sizes = genome_files.sizes
        args.pb_cfg = pipeline_args.pb_cfg
        args.chunksize = pipeline_args.pb_genome_chunksize
        args.overlap = pipeline_args.pb_genome_overlap
        args.species = pipeline_args.augustus_species
        args.hints_gff = BuildDb.get_args(pipeline_args, genome).hints_path
        args.augustus_pb_gtf = os.path.join(base_dir, genome + '.augPB.gtf')
        args.augustus_pb_gp = os.path.join(base_dir, genome + '.augPB.gp')
        args.augustus_pb_raw_gtf = os.path.join(base_dir, genome + '.raw.augPB.gtf')
        # invert the UTR flag
        args.utr = not pipeline_args.augustus_utr_off
        return args

    def validate(self):
        for tool in ['augustus', 'joingenes']:
            if not tools.misc.is_exec(tool):
                raise ToolMissingException('Auxiliary program {} from the Augustus package not in path.'.format(tool))

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.isoseq_genomes:
            yield self.clone(AugustusPbDriverTask, genome=target_genome)


class AugustusPbDriverTask(ToilTask):
    """
    Task for per-genome launching of a toil pipeline for running AugustusPB.
    """
    genome = luigi.Parameter()

    def output(self):
        pipeline_args = self.get_pipeline_args()
        augustus_pb_args = AugustusPb.get_args(pipeline_args, self.genome)
        yield luigi.LocalTarget(augustus_pb_args.augustus_pb_gp)
        yield luigi.LocalTarget(augustus_pb_args.augustus_pb_gtf)
        yield luigi.LocalTarget(augustus_pb_args.augustus_pb_raw_gtf)

    def requires(self):
        return self.clone(TransMap), self.clone(BuildDb)

    def run(self):
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'augustus_pb', self.genome)
        logger.info('Launching AugustusPB toil pipeline on {}.'.format(self.genome))
        toil_options = self.prepare_toil_options(toil_work_dir)
        augustus_pb_args = self.get_module_args(AugustusPb, genome=self.genome)
        augustus_pb(augustus_pb_args, toil_options)
        if 'stats_path' in augustus_pb_args:
            self.get_stats(toil_options, augustus_pb_args.stat_file)
        logger.info('Finished AugustusPB toil pipeline on {}.'.format(self.genome))


class FindDenovoParents(PipelineTask):
    """Task for finding parental gene candidates for denovo predictions. Flags possible fusions"""
    mode = luigi.Parameter()

    @staticmethod
    def get_args(pipeline_args, mode):
        args = tools.misc.HashableNamespace()
        if mode == 'augPB':
            args.tablename = tools.sqlInterface.AugPbAlternativeGenes.__tablename__
            args.gps = {genome: AugustusPb.get_args(pipeline_args, genome).augustus_pb_gp
                        for genome in pipeline_args.isoseq_genomes}
            args.filtered_tm_gps = {genome: TransMap.get_args(pipeline_args, genome).filtered_tm_gp
                                    for genome in pipeline_args.isoseq_genomes - {pipeline_args.ref_genome}}
            args.unfiltered_tm_gps = {genome: TransMap.get_args(pipeline_args, genome).tm_gp
                                      for genome in pipeline_args.isoseq_genomes - {pipeline_args.ref_genome}}
            args.chrom_sizes = {genome: GenomeFiles.get_args(pipeline_args, genome).sizes
                                for genome in pipeline_args.isoseq_genomes}
            # add the reference annotation as a pseudo-transMap to assign parents in reference
            args.filtered_tm_gps[pipeline_args.ref_genome] = ReferenceFiles.get_args(pipeline_args).annotation_gp
            args.unfiltered_tm_gps[pipeline_args.ref_genome] = ReferenceFiles.get_args(pipeline_args).annotation_gp
        elif mode == 'augCGP':
            args.tablename = tools.sqlInterface.AugCgpAlternativeGenes.__tablename__
            args.gps = AugustusCgp.get_args(pipeline_args).augustus_cgp_gp
            filtered_tm_gp_files = {genome: TransMap.get_args(pipeline_args, genome).filtered_tm_gp
                                    for genome in pipeline_args.target_genomes}
            unfiltered_tm_gp_files = {genome: TransMap.get_args(pipeline_args, genome).tm_gp
                                      for genome in pipeline_args.target_genomes}
            # add the reference annotation as a pseudo-transMap to assign parents in reference
            filtered_tm_gp_files[pipeline_args.ref_genome] = ReferenceFiles.get_args(pipeline_args).annotation_gp
            unfiltered_tm_gp_files[pipeline_args.ref_genome] = ReferenceFiles.get_args(pipeline_args).annotation_gp
            args.filtered_tm_gps = filtered_tm_gp_files
            args.unfiltered_tm_gps = unfiltered_tm_gp_files
            args.chrom_sizes = {genome: GenomeFiles.get_args(pipeline_args, genome).sizes
                                for genome in list(pipeline_args.target_genomes) + [pipeline_args.ref_genome]}
        else:
            raise Exception('Invalid mode passed to FindDenovoParents')
        return args

    def requires(self):
        if self.mode == 'augPB':
            yield self.clone(AugustusPb)
        elif self.mode == 'augCGP':
            yield self.clone(AugustusCgp)
        else:
            raise Exception('Invalid mode passed to FindDenovoParents')
        yield self.clone(TransMap)

    def get_table_targets(self, genome, tablename, pipeline_args):
        db = pipeline_args.dbs[genome]
        tools.fileOps.ensure_file_dir(db)
        conn_str = 'sqlite:///{}'.format(db)
        return luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                   target_table=tablename,
                                                   update_id='_'.join([tablename, str(hash(pipeline_args))]))

    def output(self):
        pipeline_args = self.get_pipeline_args()
        denovo_args = FindDenovoParents.get_args(pipeline_args, self.mode)
        for genome in denovo_args.gps:
            yield self.get_table_targets(genome, denovo_args.tablename, pipeline_args)

    def run(self):
        pipeline_args = self.get_pipeline_args()
        denovo_args = FindDenovoParents.get_args(pipeline_args, self.mode)
        for genome, denovo_gp in denovo_args.gps.iteritems():
            table_target = self.get_table_targets(genome, denovo_args.tablename, pipeline_args)
            filtered_tm_gp = denovo_args.filtered_tm_gps[genome]
            unfiltered_tm_gp = denovo_args.unfiltered_tm_gps[genome]
            chrom_sizes = denovo_args.chrom_sizes[genome]
            df = assign_parents(filtered_tm_gp, unfiltered_tm_gp, chrom_sizes, denovo_gp)
            db = pipeline_args.dbs[genome]
            with tools.sqlite.ExclusiveSqlConnection(db) as engine:
                df.to_sql(denovo_args.tablename, engine, if_exists='replace')
            table_target.touch()
            counts = collections.Counter(df.ResolutionMethod)
            log_msg = 'Loaded table: {}.{}. Results: {}'
            assigned_str = '{}: {:,}'.format('assigned', counts[None])
            log_msg = log_msg.format(genome, denovo_args.tablename, assigned_str)
            result_str = ', '.join(['{}: {:,}'.format(name, val)
                                    for name, val in sorted(counts.iteritems()) if name is not None])
            if len(result_str) > 0:
                log_msg += ', ' + result_str + '.'
            logger.info(log_msg)


class Hgm(PipelineWrapperTask):
    """
    Task for launching the HomGeneMapping toil pipeline. This pipeline finds the cross species RNA-seq and annotation
    support across all species.
    It will be launched once for each of transMap, AugustusTM, AugustusTMR, AugustusCGP
    """
    @staticmethod
    def get_args(pipeline_args, mode):
        base_dir = os.path.join(pipeline_args.work_dir, 'hgm', mode)
        if mode == 'augCGP':
            # add reference to the target genomes
            tgt_genomes = list(pipeline_args.target_genomes) + [pipeline_args.ref_genome]
            gtf_in_files = {genome: AugustusCgp.get_args(pipeline_args).augustus_cgp_gtf[genome]
                            for genome in tgt_genomes}
        elif mode == 'augTM':
            tgt_genomes = pipeline_args.target_genomes
            gtf_in_files = {genome: Augustus.get_args(pipeline_args, genome).augustus_tm_gtf
                            for genome in tgt_genomes}
        elif mode == 'augTMR':
            # remove reference it may have RNA-seq
            tgt_genomes = (pipeline_args.rnaseq_genomes & set(pipeline_args.target_genomes)) - {pipeline_args.ref_genome}
            gtf_in_files = {genome: Augustus.get_args(pipeline_args, genome).augustus_tmr_gtf
                            for genome in tgt_genomes}
        elif mode == 'augPB':
            # add reference genome to target_genomes, but then intersect with isoseq genomes
            tgt_genomes = (set(pipeline_args.target_genomes) | {pipeline_args.ref_genome}) & pipeline_args.isoseq_genomes
            gtf_in_files = {genome: AugustusPb.get_args(pipeline_args, genome).augustus_pb_gtf
                            for genome in tgt_genomes}
        elif mode == 'transMap':
            tgt_genomes = pipeline_args.target_genomes
            gtf_in_files = {genome: TransMap.get_args(pipeline_args, genome).tm_gtf
                            for genome in tgt_genomes}
        else:
            raise UserException('Invalid mode was passed to Hgm module: {}.'.format(mode))
        args = tools.misc.HashableNamespace()
        args.genomes = tgt_genomes
        args.ref_genome = pipeline_args.ref_genome
        args.hal = pipeline_args.hal
        args.in_gtf = gtf_in_files
        args.gtf_out_dir = base_dir
        args.gtf_out_files = {genome: os.path.join(base_dir, genome + '.gtf') for genome in tgt_genomes}
        args.hints_db = pipeline_args.hints_db
        args.annotation_gtf = ReferenceFiles.get_args(pipeline_args).annotation_gtf
        args.annotation_gp = ReferenceFiles.get_args(pipeline_args).annotation_gp
        args.hgm_cpu = pipeline_args.hgm_cpu
        return args

    def validate(self):
        for tool in ['homGeneMapping', 'join_mult_hints.pl']:
            if not tools.misc.is_exec(tool):
                raise ToolMissingException('auxiliary program {} from the Augustus '
                                           'package not in global path.'.format(tool))
        if not tools.misc.is_exec('halLiftover'):
            raise ToolMissingException('halLiftover from the halTools package not in global path.')
        for tool in ['bedtools', 'bedSort']:
            if not tools.misc.is_exec(tool):
                raise ToolMissingException('{} is required for the homGeneMapping module.'.format(tool))

    def requires(self):
        pipeline_args = self.get_pipeline_args()
        self.validate()
        for mode in pipeline_args.modes:
            yield self.clone(HgmDriverTask, mode=mode)


class HgmDriverTask(PipelineTask):
    """
    Task for running each individual instance of the Hgm pipeline. Dumps the results into a sqlite database
    Also produces a GTF file that is parsed into this database.
    """
    mode = luigi.Parameter()

    def output(self):
        pipeline_args = self.get_pipeline_args()
        hgm_args = Hgm.get_args(pipeline_args, self.mode)
        for genome in hgm_args.genomes:
            db = pipeline_args.dbs[genome]
            tools.fileOps.ensure_file_dir(db)
            conn_str = 'sqlite:///{}'.format(db)
            tablename = tools.sqlInterface.tables['hgm'][self.mode].__tablename__
            yield luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                      target_table=tablename,
                                                      update_id='_'.join([tablename, str(hash(pipeline_args))]))
        for f in hgm_args.gtf_out_files.itervalues():
            yield luigi.LocalTarget(f)

    def requires(self):
        if self.mode == 'augCGP':
            yield self.clone(AugustusCgp)
            yield self.clone(FindDenovoParents, mode='augCGP')
        elif self.mode == 'augTM' or self.mode == 'augTMR':
            yield self.clone(Augustus)
        elif self.mode == 'transMap':
            yield self.clone(TransMap)
        elif self.mode == 'augPB':
            yield self.clone(AugustusPb)
            yield self.clone(FindDenovoParents, mode='augPB')
        else:
            raise UserException('Invalid mode passed to HgmDriverTask: {}.'.format(self.mode))
        yield self.clone(BuildDb)
        yield self.clone(ReferenceFiles)

    def run(self):
        logger.info('Launching homGeneMapping for {}.'.format(self.mode))
        pipeline_args = self.get_pipeline_args()
        hgm_args = Hgm.get_args(pipeline_args, self.mode)
        hgm(hgm_args)
        # convert the output to a dataframe and write to the genome database
        databases = self.__class__.get_databases(pipeline_args)
        tablename = tools.sqlInterface.tables['hgm'][self.mode].__tablename__
        for genome, sqla_target in itertools.izip(*[hgm_args.genomes, self.output()]):
            df = parse_hgm_gtf(hgm_args.gtf_out_files[genome], genome)
            with tools.sqlite.ExclusiveSqlConnection(databases[genome]) as engine:
                df.to_sql(tablename, engine, if_exists='replace')
            sqla_target.touch()
            logger.info('Loaded table: {}.{}'.format(genome, tablename))


class IsoSeqTranscripts(PipelineWrapperTask):
    """
    Constructs a database table representing all unique exon structures seen in an IsoSeq dataset. This is used
    to validate isoforms (compared to homGeneMapping, which just validates exons/introns independently).

    These structures are analogous to Transcript objects.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        args = tools.misc.HashableNamespace()
        args.genome = genome
        args.hints_gff = BuildDb.get_args(pipeline_args, genome).hints_path
        return args

    def requires(self):
        pipeline_args = self.get_pipeline_args()
        for genome in pipeline_args.isoseq_genomes:
            yield self.clone(IsoSeqTranscriptsDriverTask, genome=genome)


class IsoSeqTranscriptsDriverTask(PipelineTask):
    """
    Driver task for IsoSeqTranscripts
    """
    genome = luigi.Parameter()
    tablename = tools.sqlInterface.IsoSeqExonStructures.__tablename__

    def output(self):
        pipeline_args = self.get_pipeline_args()
        db = pipeline_args.dbs[self.genome]
        tools.fileOps.ensure_file_dir(db)
        conn_str = 'sqlite:///{}'.format(db)
        return luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                   target_table=self.tablename,
                                                   update_id='_'.join([self.tablename, str(hash(pipeline_args))]))

    def requires(self):
        yield self.clone(BuildDb)

    def construct_intervals(self, hints):
        """
        Converts hints derived from IsoSeq BAMs into discrete clusters of transcript objects. Merges all alignment gaps
        below 50bp and separates clusters over 100kb separated to avoid mega-transcripts for tandem gene families.
        """
        lines = [x.split() for x in open(hints) if 'PB' in x and '\texon\t' in x]
        # group these exons by grp tag
        groups = collections.defaultdict(list)
        for l in lines:
            attrs = dict([x.split('=') for x in l[-1].split(';')])
            if 'grp' not in attrs:  # not all introns get confidently assigned a group
                continue
            groups[attrs['grp']].append([l[0], int(l[3]) - 1, int(l[4])])

        # for each grp, perform clustering with 100kb distance to separate contigs as well as disjoint mappings
        # to do this, we use the ClusterTree data structure from bx-python.
        # we need one cluster-tree for every chromosome-grp combination
        # the 1 after the 100000 says 'read depth of 1 is sufficient for output'
        # see https://bcbio.wordpress.com/2009/04/29/finding-and-displaying-short-reads-clustered-in-the-genome/
        cluster_trees = collections.defaultdict(lambda: collections.defaultdict(lambda: ClusterTree(100000, 1)))
        # we also need to keep track of every interval for downstream processing
        i = 0
        interval_flat_list = []
        for grp, intervals in groups.iteritems():
            for chrom, start, stop in intervals:
                interval_flat_list.append([chrom, start, stop])
                cluster_trees[chrom][grp].insert(start, stop, i)
                i += 1

        # for each cluster, convert to a transcript object
        txs = []
        for chrom in cluster_trees:
            for grp, cluster_tree in cluster_trees[chrom].iteritems():
                for start, end, interval_indices in cluster_tree.getregions():
                    intervals = [interval_flat_list[i] for i in interval_indices]
                    intervals = {tools.intervals.ChromosomeInterval(chrom, start, stop, '.')
                                 for chrom, start, stop in intervals}
                    intervals = tools.intervals.gap_merge_intervals(intervals, 50)
                    txs.append(tools.transcripts.intervals_to_bed(intervals, name=grp))

        # convert these to a dataframe for sql output
        txs = [x.get_bed() for x in txs]
        df = pd.DataFrame(txs, columns=['chromosome', 'start', 'stop', 'name', 'score', 'strand', 'thickStart',
                                        'thickStop', 'rgb', 'blockCount', 'blockSizes', 'blockStarts'])
        return df

    def run(self):
        pipeline_args = self.get_pipeline_args()
        intron_args = IsoSeqTranscripts.get_args(pipeline_args, self.genome)
        df = pd.DataFrame(self.construct_intervals(intron_args.hints_gff))
        with tools.sqlite.ExclusiveSqlConnection(pipeline_args.dbs[self.genome]) as engine:
            df.to_sql(self.tablename, engine, if_exists='replace')
        self.output().touch()
        logger.info('Loaded table {}.{}'.format(self.genome, self.tablename))


class AlignTranscripts(PipelineWrapperTask):
    """
    Aligns the transcripts from transMap/AugustusTMR to the parent transcript(s).
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'transcript_alignment')
        args = tools.misc.HashableNamespace()
        args.ref_genome = pipeline_args.ref_genome
        args.genome = genome
        args.ref_genome_fasta = GenomeFiles.get_args(pipeline_args, pipeline_args.ref_genome).fasta
        args.genome_fasta = GenomeFiles.get_args(pipeline_args, genome).fasta
        args.annotation_gp = ReferenceFiles.get_args(pipeline_args).annotation_gp
        args.ref_db_path = PipelineTask.get_database(pipeline_args, pipeline_args.ref_genome)
        # the alignment_modes members hold the input genePreds and the mRNA/CDS alignment output paths
        args.transcript_modes = {'transMap': {'gp': TransMap.get_args(pipeline_args, genome).filtered_tm_gp,
                                              'mRNA': os.path.join(base_dir, genome + '.transMap.mRNA.psl'),
                                              'CDS': os.path.join(base_dir, genome + '.transMap.CDS.psl')}}
        if pipeline_args.augustus is True:
            args.transcript_modes['augTM'] = {'gp':  Augustus.get_args(pipeline_args, genome).augustus_tm_gp,
                                              'mRNA': os.path.join(base_dir, genome + '.augTM.mRNA.psl'),
                                              'CDS': os.path.join(base_dir, genome + '.augTM.CDS.psl')}
        if pipeline_args.augustus is True and genome in pipeline_args.rnaseq_genomes:
            args.transcript_modes['augTMR'] = {'gp': Augustus.get_args(pipeline_args, genome).augustus_tmr_gp,
                                               'mRNA': os.path.join(base_dir, genome + '.augTMR.mRNA.psl'),
                                               'CDS': os.path.join(base_dir, genome + '.augTMR.CDS.psl')}
        return args

    def validate(self):
        for tool in ['blat', 'pslCheck']:
            if not tools.misc.is_exec(tool):
                raise ToolMissingException('Tool {} not in global path.'.format(tool))

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            yield self.clone(AlignTranscriptDriverTask, genome=target_genome)


class AlignTranscriptDriverTask(ToilTask):
    """
    Task for per-genome launching of a toil pipeline for aligning all transcripts found back to the reference in
    transcript space using BLAT.

    Each task returns a PSL of all alignments that will be analyzed next by EvaluateTranscripts.
    """
    genome = luigi.Parameter()

    def output(self):
        alignment_args = self.get_module_args(AlignTranscripts, genome=self.genome)
        for mode, paths in alignment_args.transcript_modes.iteritems():
            for aln_type in ['CDS', 'mRNA']:
                yield luigi.LocalTarget(paths[aln_type])

    def requires(self):
        alignment_args = self.get_module_args(AlignTranscripts, genome=self.genome)
        if 'augTM' in alignment_args.transcript_modes:
            yield self.clone(Augustus)
        yield self.clone(TransMap)
        yield self.clone(ReferenceFiles)
        yield self.clone(GenomeFiles)

    def run(self):
        logger.info('Launching Align Transcript toil pipeline for {} using {}.'.format(self.genome, self.batchSystem))
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'transcript_alignment', self.genome)
        toil_options = self.prepare_toil_options(toil_work_dir)
        alignment_args = self.get_module_args(AlignTranscripts, genome=self.genome)
        align_transcripts(alignment_args, toil_options)
        logger.info('Align Transcript toil pipeline for {} completed.'.format(self.genome))


class EvaluateTranscripts(PipelineWrapperTask):
    """
    Evaluates all transcripts for important features. See the classify.py module for details on how this works.

    Each task will generate a genome-specific sqlite database. See the classify.py docstring for details.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        args = tools.misc.HashableNamespace()
        args.db_path = pipeline_args.dbs[genome]
        args.ref_db_path = PipelineTask.get_database(pipeline_args, pipeline_args.ref_genome)
        args.annotation_gp = ReferenceFiles.get_args(pipeline_args).annotation_gp
        args.fasta = GenomeFiles.get_args(pipeline_args, genome).fasta
        args.genome = genome
        args.ref_genome = pipeline_args.ref_genome
        # pass along all of the paths from alignment
        args.transcript_modes = AlignTranscripts.get_args(pipeline_args, genome).transcript_modes
        return args

    def validate(self):
        pass

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            yield self.clone(EvaluateDriverTask, genome=target_genome)


class EvaluateDriverTask(PipelineTask):
    """
    Task for per-genome launching of a toil pipeline for aligning transcripts to their parent.
    """
    genome = luigi.Parameter()

    def build_table_names(self, eval_args):
        """construct table names based on input arguments"""
        tables = []
        for aln_mode in ['mRNA', 'CDS']:
            for tx_mode in eval_args.transcript_modes.iterkeys():
                names = [x.__tablename__ for x in tools.sqlInterface.tables[aln_mode][tx_mode].values()]
                tables.extend(names)
        return tables

    def pair_table_output(self, eval_args):
        """return dict of {table_name: SQLAlchemyTarget} for final writing"""
        return dict(zip(*[self.build_table_names(eval_args), self.output()]))

    def write_to_sql(self, results, eval_args):
        """Load the results into the SQLite database"""
        with tools.sqlite.ExclusiveSqlConnection(eval_args.db_path) as engine:
            for table, target in self.pair_table_output(eval_args).iteritems():
                df = results[table]
                df.to_sql(table, engine, if_exists='replace')
                target.touch()
                logger.info('Loaded table: {}.{}'.format(self.genome, table))

    def output(self):
        pipeline_args = self.get_pipeline_args()
        eval_args = self.get_module_args(EvaluateTranscripts, genome=self.genome)
        tools.fileOps.ensure_file_dir(eval_args.db_path)
        conn_str = 'sqlite:///{}'.format(eval_args.db_path)
        for table in self.build_table_names(eval_args):
            yield luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                      target_table=table,
                                                      update_id='_'.join([table, str(hash(pipeline_args))]))

    def requires(self):
        return self.clone(AlignTranscripts), self.clone(ReferenceFiles), self.clone(TransMap)

    def run(self):
        logger.info('Evaluating transcript alignments for {}.'.format(self.genome))
        eval_args = self.get_module_args(EvaluateTranscripts, genome=self.genome)
        results = classify(eval_args)
        # results should be a dictionary of {table: dataframe}
        self.write_to_sql(results, eval_args)


class Consensus(PipelineWrapperTask):
    """
    Construct the consensus gene sets making use of the classification databases.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.out_dir, 'consensus_gene_set')
        # grab the genePred of every mode
        args = tools.misc.HashableNamespace()
        gp_list = [TransMap.get_args(pipeline_args, genome).filtered_tm_gp]
        args.tx_modes = ['transMap']
        args.denovo_tx_modes = []
        if pipeline_args.augustus is True:
            gp_list.append(Augustus.get_args(pipeline_args, genome).augustus_tm_gp)
            args.tx_modes.append('augTM')
        if pipeline_args.augustus is True and genome in pipeline_args.rnaseq_genomes:
            gp_list.append(Augustus.get_args(pipeline_args, genome).augustus_tmr_gp)
            args.tx_modes.append('augTMR')
        if pipeline_args.augustus_cgp is True:
            gp_list.append(AugustusCgp.get_args(pipeline_args).augustus_cgp_gp[genome])
            args.denovo_tx_modes.append('augCGP')
        if pipeline_args.augustus_pb is True and genome in pipeline_args.isoseq_genomes:
            gp_list.append(AugustusPb.get_args(pipeline_args, genome).augustus_pb_gp)
            args.denovo_tx_modes.append('augPB')
        args.gp_list = gp_list
        args.genome = genome
        args.transcript_modes = AlignTranscripts.get_args(pipeline_args, genome).transcript_modes.keys()
        args.augustus_cgp = pipeline_args.augustus_cgp
        args.db_path = pipeline_args.dbs[genome]
        args.ref_db_path = PipelineTask.get_database(pipeline_args, pipeline_args.ref_genome)
        args.hints_db_has_rnaseq = len(pipeline_args.rnaseq_genomes) > 0
        args.annotation_gp = ReferenceFiles.get_args(pipeline_args).annotation_gp
        args.fasta = GenomeFiles.get_args(pipeline_args, genome).fasta
        args.consensus_gp = os.path.join(base_dir, genome + '.gp')
        args.consensus_gp_info = os.path.join(base_dir, genome + '.gp_info')
        args.consensus_gff3 = os.path.join(base_dir, genome + '.gff3')
        args.consensus_fasta = os.path.join(base_dir, genome + '.consensus.fasta')
        args.consensus_protein_fasta = os.path.join(base_dir, genome + '.protein.consensus.fasta')
        args.metrics_json = os.path.join(PipelineTask.get_metrics_dir(pipeline_args, genome), 'consensus.json')
        # user configurable options on how consensus finding should work
        args.intron_rnaseq_support = pipeline_args.intron_rnaseq_support
        args.exon_rnaseq_support = pipeline_args.exon_rnaseq_support
        args.intron_annot_support = pipeline_args.intron_annot_support
        args.exon_annot_support = pipeline_args.exon_annot_support
        args.original_intron_support = pipeline_args.original_intron_support
        args.denovo_num_introns = pipeline_args.denovo_num_introns
        args.denovo_splice_support = pipeline_args.denovo_splice_support
        args.denovo_exon_support = pipeline_args.denovo_exon_support
        args.require_pacbio_support = pipeline_args.require_pacbio_support
        args.in_species_rna_support_only = pipeline_args.in_species_rna_support_only
        return args

    def validate(self):
        pass

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            yield self.clone(ConsensusDriverTask, genome=target_genome)


class ConsensusDriverTask(RebuildableTask):
    """
    Driver task for performing consensus finding.
    """
    genome = luigi.Parameter()

    def output(self):
        consensus_args = self.get_module_args(Consensus, genome=self.genome)
        yield luigi.LocalTarget(consensus_args.metrics_json)
        yield luigi.LocalTarget(consensus_args.consensus_gp)
        yield luigi.LocalTarget(consensus_args.consensus_gp_info)
        yield luigi.LocalTarget(consensus_args.consensus_gff3)
        yield luigi.LocalTarget(consensus_args.consensus_fasta)
        yield luigi.LocalTarget(consensus_args.consensus_protein_fasta)

    def requires(self):
        pipeline_args = self.get_pipeline_args()
        yield self.clone(EvaluateTransMap)
        yield self.clone(EvaluateTranscripts)
        yield self.clone(Hgm)
        yield self.clone(AlignTranscripts)
        yield self.clone(ReferenceFiles)
        yield self.clone(EvaluateTransMap)
        yield self.clone(TransMap)
        if pipeline_args.augustus_pb:
            yield self.clone(AugustusPb)
            yield self.clone(IsoSeqTranscripts)
            yield self.clone(FindDenovoParents, mode='augPB')
        if pipeline_args.augustus_cgp:
            yield self.clone(AugustusCgp)
            yield self.clone(FindDenovoParents, mode='augCGP')

    def run(self):
        consensus_args = self.get_module_args(Consensus, genome=self.genome)
        logger.info('Generating consensus gene set for {}.'.format(self.genome))
        metrics_dict = generate_consensus(consensus_args)
        metrics_json = self.output().next()
        PipelineTask.write_metrics(metrics_dict, metrics_json)


class Plots(RebuildableTask):
    """
    Produce final analysis plots
    """
    @staticmethod
    def get_args(pipeline_args):
        base_dir = os.path.join(pipeline_args.out_dir, 'plots')
        ordered_genomes = tools.hal.build_genome_order(pipeline_args.hal, pipeline_args.ref_genome,
                                                       pipeline_args.target_genomes,
                                                       pipeline_args.annotate_ancestors)
        args = tools.misc.HashableNamespace()
        args.ordered_genomes = ordered_genomes
        # plots derived from transMap results
        args.tm_coverage = luigi.LocalTarget(os.path.join(base_dir, 'transmap_coverage.pdf'))
        args.tm_identity = luigi.LocalTarget(os.path.join(base_dir, 'transmap_identity.pdf'))
        # plots derived from transMap filtering
        args.paralogy = luigi.LocalTarget(os.path.join(base_dir, 'paralogy.pdf'))
        args.gene_collapse = luigi.LocalTarget(os.path.join(base_dir, 'gene_family_collapse.pdf'))
        # plots derived from transcript alignment / consensus finding
        args.coverage = luigi.LocalTarget(os.path.join(base_dir, 'coverage.pdf'))
        args.identity = luigi.LocalTarget(os.path.join(base_dir, 'identity.pdf'))
        args.completeness = luigi.LocalTarget(os.path.join(base_dir, 'completeness.pdf'))
        args.consensus_extrinsic_support = luigi.LocalTarget(os.path.join(base_dir, 'consensus_extrinsic_support.pdf'))
        args.consensus_annot_support = luigi.LocalTarget(os.path.join(base_dir, 'consensus_annotation_support.pdf'))
        args.tx_modes = luigi.LocalTarget(os.path.join(base_dir, 'transcript_modes.pdf'))
        args.indel = luigi.LocalTarget(os.path.join(base_dir, 'coding_indels.pdf'))
        args.missing = luigi.LocalTarget(os.path.join(base_dir, 'missing_genes_transcripts.pdf'))
        # plots that depend on execution mode
        if pipeline_args.augustus is True:
            args.improvement = luigi.LocalTarget(os.path.join(base_dir, 'augustus_improvement.pdf'))
        if 'augCGP' in pipeline_args.modes or 'augPB' in pipeline_args.modes:
            args.denovo = luigi.LocalTarget(os.path.join(base_dir, 'denovo.pdf'))
        if 'augPB' in pipeline_args.modes:
            args.pb_support = luigi.LocalTarget(os.path.join(base_dir, 'IsoSeq_isoform_validation.pdf'))
            args.pb_genomes = pipeline_args.isoseq_genomes
        args.split_genes = luigi.LocalTarget(os.path.join(base_dir, 'split_genes.pdf'))
        # input data
        args.metrics_jsons = OrderedDict([[genome, Consensus.get_args(pipeline_args, genome).metrics_json]
                                          for genome in ordered_genomes])
        args.tm_jsons = OrderedDict([[genome, TransMap.get_args(pipeline_args, genome).metrics_json]
                                     for genome in ordered_genomes])
        args.annotation_db = PipelineTask.get_database(pipeline_args, pipeline_args.ref_genome)
        args.dbs = OrderedDict([[genome, PipelineTask.get_database(pipeline_args, genome)]
                                for genome in ordered_genomes])
        args.in_species_rna_support_only = pipeline_args.in_species_rna_support_only
        return args

    def output(self):
        pipeline_args = self.get_pipeline_args()
        args = Plots.get_args(pipeline_args)
        return [p for p in args.__dict__.itervalues() if isinstance(p, luigi.LocalTarget)]

    def requires(self):
        yield self.clone(Consensus)
        yield self.clone(EvaluateTransMap)
        yield self.clone(EvaluateTranscripts)
        yield self.clone(Hgm)
        yield self.clone(ReferenceFiles)
        yield self.clone(EvaluateTransMap)
        yield self.clone(TransMap)

    def run(self):
        pipeline_args = self.get_pipeline_args()
        logger.info('Generating plots.')
        generate_plots(Plots.get_args(pipeline_args))


class ReportStats(PipelineTask):
    """
    Reports all the stats at the end of the pipeline
    """
    def requires(self):
        yield self.clone(PrepareFiles)
        yield self.clone(BuildDb)
        yield self.clone(Chaining)
        yield self.clone(TransMap)
        yield self.clone(EvaluateTransMap)
        if self.augustus is True:
            yield self.clone(Augustus)
        if self.augustus_cgp is True:
            yield self.clone(AugustusCgp)
            yield self.clone(FindDenovoParents, mode='augCGP')
        if self.augustus_pb is True:
            yield self.clone(AugustusPb)
            yield self.clone(FindDenovoParents, mode='augPB')
            yield self.clone(IsoSeqTranscripts)
        yield self.clone(Hgm)
        yield self.clone(AlignTranscripts)
        yield self.clone(EvaluateTranscripts)
        yield self.clone(Consensus)
        yield self.clone(Plots)
        if self.assembly_hub is True:
            yield self.clone(AssemblyHub)

    def output(self):
        # dumb -- need it to be something
        pipeline_args = self.get_pipeline_args()
        tools.fileOps.ensure_file_dir(pipeline_args.stats_db)
        conn_str = 'sqlite:///{}'.format(pipeline_args.stats_db)
        return luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                   target_table='stats',
                                                   update_id='_'.join(['stats', str(hash(pipeline_args))]))

    def run(self):
        pipeline_args = self.get_pipeline_args()
        luigi_stats = tools.sqlInterface.load_luigi_stats(pipeline_args.stats_db, 'stats')
        
        try:
            toil_stats = tools.sqlInterface.load_luigi_stats(pipeline_args.stats_db, 'toil_stats')
        except ValueError:
            logger.warning('Toil task already ran, therefore no stats')
        else:
            core_time = round(sum(luigi_stats.ProcessingTime) / 3600, 1)
            toil_core_time = round(sum(toil_stats.TotalTime) / 3600, 1)
            total = core_time + toil_core_time
            logger.info('Local core time: {:,} hours. Toil core time: {:,} hours. '
                    'Total computation time: {:,} hours.'.format(core_time, toil_core_time, total))
        self.output().touch()


class AssemblyHub(PipelineWrapperTask):
    """
    Construct an assembly hub out of all the results
    """
    def requires(self):
        tools.fileOps.ensure_dir(self.out_dir)
        yield self.clone(CreateDirectoryStructure)
        yield self.clone(CreateTracks)


class CreateDirectoryStructure(RebuildableTask):
    """
    Constructs the directory structure. Creates symlinks for all relevant files.
    """
    @staticmethod
    def get_args(pipeline_args):
        args = tools.misc.HashableNamespace()
        args.genomes = list(pipeline_args.target_genomes) + [pipeline_args.ref_genome]
        args.out_dir = os.path.join(pipeline_args.out_dir, 'assemblyHub')
        args.hub_txt = os.path.join(args.out_dir, 'hub.txt')
        args.genomes_txt = os.path.join(args.out_dir, 'genomes.txt')
        args.groups_txt = os.path.join(args.out_dir, 'groups.txt')
        genome_files = frozendict({genome: GenomeFiles.get_args(pipeline_args, genome) for genome in args.genomes})
        sizes = {}
        twobits = {}
        trackdbs = {}
        for genome, genome_file in genome_files.iteritems():
            sizes[genome] = (genome_file.sizes, os.path.join(args.out_dir, genome, 'chrom.sizes'))
            twobits[genome] = (genome_file.two_bit, os.path.join(args.out_dir, genome, '{}.2bit'.format(genome)))
            trackdbs[genome] = os.path.join(args.out_dir, genome, 'trackDb.txt')
        args.sizes = frozendict(sizes)
        args.twobits = frozendict(twobits)
        args.trackdbs = frozendict(trackdbs)
        args.hal = os.path.join(args.out_dir, os.path.basename(pipeline_args.hal))
        return args

    def requires(self):
        yield self.clone(Consensus)
        yield self.clone(EvaluateTransMap)
        yield self.clone(EvaluateTranscripts)
        yield self.clone(Hgm)
        yield self.clone(ReferenceFiles)
        yield self.clone(EvaluateTransMap)
        yield self.clone(TransMap)

    def output(self):
        pipeline_args = self.get_pipeline_args()
        args = CreateDirectoryStructure.get_args(pipeline_args)
        yield luigi.LocalTarget(args.hub_txt)
        yield luigi.LocalTarget(args.genomes_txt)
        yield luigi.LocalTarget(args.groups_txt)
        yield luigi.LocalTarget(args.hal)
        for local_path, hub_path in args.sizes.itervalues():
            yield luigi.LocalTarget(hub_path)
        for local_path, hub_path in args.twobits.itervalues():
            yield luigi.LocalTarget(hub_path)

    def run(self):
        pipeline_args = self.get_pipeline_args()
        args = CreateDirectoryStructure.get_args(pipeline_args)
        tools.fileOps.ensure_file_dir(args.out_dir)

        # write the hub.txt file
        with luigi.LocalTarget(args.hub_txt).open('w') as outf:
            outf.write(hub_str.format(hal=os.path.splitext(os.path.basename(pipeline_args.hal))[0],
                                      email=pipeline_args.hub_email))

        # write the groups.txt file
        with luigi.LocalTarget(args.groups_txt).open('w') as outf:
            outf.write(groups_str)

        # write the genomes.txt file, construct a dir
        with luigi.LocalTarget(args.genomes_txt).open('w') as outf:
            for genome, (sizes_local_path, sizes_hub_path) in args.sizes.iteritems():
                outf.write(genome_str.format(genome=genome, default_pos=find_default_pos(sizes_local_path)))

        # link the hal
        shutil.copy(pipeline_args.hal, args.hal)

        # construct a directory for each genome
        for genome, (sizes_local_path, sizes_hub_path) in args.sizes.iteritems():
            tools.fileOps.ensure_file_dir(sizes_hub_path)
            shutil.copy(sizes_local_path, sizes_hub_path)
            twobit_local_path, twobit_hub_path = args.twobits[genome]
            shutil.copy(twobit_local_path, twobit_hub_path)


class CreateTracks(PipelineWrapperTask):
    """
    Wrapper task for track creation.
    """
    def validate(self):
        for tool in ['bedSort', 'pslToBigPsl', 'wiggletools', 'wigToBigWig']:#, 'bamCoverage']:
            if not tools.misc.is_exec(tool):
                raise ToolMissingException('Tool {} not in global path.'.format(tool))

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for genome in pipeline_args.target_genomes:
            yield self.clone(CreateTracksDriverTask, genome=genome)
            yield self.clone(CreateTrackDbs, genome=genome)
        yield self.clone(CreateTracksDriverTask, genome=pipeline_args.ref_genome)
        yield self.clone(CreateTrackDbs, genome=pipeline_args.ref_genome)


class CreateTracksDriverTask(PipelineWrapperTask):
    """
    Dynamically generates each track task, the combines the results into the final trackDb.
    """
    genome = luigi.Parameter()

    def requires(self):
        pipeline_args = self.get_pipeline_args()
        if self.genome not in pipeline_args.target_genomes and self.genome != pipeline_args.ref_genome:
            return
        directory_args = CreateDirectoryStructure.get_args(pipeline_args)
        out_dir = os.path.join(directory_args.out_dir, self.genome)
        if pipeline_args.augustus_cgp is True and self.genome in pipeline_args.target_genomes:
            yield self.clone(DenovoTrack, track_path=os.path.join(out_dir, 'augustus_cgp.bb'),
                             trackdb_path=os.path.join(out_dir, 'augustus_cgp.txt'), mode='augCGP')
        if pipeline_args.augustus_pb is True and self.genome in pipeline_args.isoseq_genomes:
            yield self.clone(DenovoTrack, track_path=os.path.join(out_dir, 'augustus_pb.bb'),
                             trackdb_path=os.path.join(out_dir, 'augustus_pb.txt'), mode='augPB')

        # TODO: allow non-reference annotations
        #if self.genome in pipeline_args.annotation_genomes:
        if self.genome == pipeline_args.ref_genome:
            annotation_gp = ReferenceFiles.get_args(pipeline_args).annotation_gp
            yield self.clone(BgpTrack, track_path=os.path.join(out_dir, 'annotation.bb'),
                             trackdb_path=os.path.join(out_dir, 'annotation.txt'),
                             genepred_path=annotation_gp, label=os.path.splitext(os.path.basename(annotation_gp))[0])

        if self.genome in pipeline_args.target_genomes:
            yield self.clone(ConsensusTrack, track_path=os.path.join(out_dir, 'consensus.bb'),
                             trackdb_path=os.path.join(out_dir, 'consensus.txt'))

            tx_modes = ['transMap']
            if pipeline_args.augustus is True:
                tx_modes.append('augTM')
                if self.genome in pipeline_args.rnaseq_genomes:
                    tx_modes.append('augTMR')
            yield self.clone(EvaluationTrack, track_path=os.path.join(out_dir, 'evaluation.bb'),
                             trackdb_path=os.path.join(out_dir, 'evaluation.txt'),
                             tx_modes=tuple(tx_modes))

            tm_args = TransMap.get_args(pipeline_args, self.genome)
            yield self.clone(TransMapTrack, track_path=os.path.join(out_dir, 'transmap.bb'),
                             trackdb_path=os.path.join(out_dir, 'transmap.txt'))
            yield self.clone(BgpTrack, track_path=os.path.join(out_dir, 'filtered_transmap.bb'),
                             trackdb_path=os.path.join(out_dir, 'filtered_transmap.txt'),
                             genepred_path=tm_args.filtered_tm_gp, label='Filtered transMap', visibility='hide')

            if pipeline_args.augustus is True and self.genome in pipeline_args.rnaseq_genomes:
                yield self.clone(AugustusTrack, track_path=os.path.join(out_dir, 'augustus.bb'),
                                 trackdb_path=os.path.join(out_dir, 'augustus.txt'))

        if self.genome in pipeline_args.isoseq_genomes:
            isoseq_bams = []
            # add a number to make names unique
            for i, bam in enumerate(pipeline_args.cfg['ISO_SEQ_BAM'][self.genome]):
                new_bam = os.path.join(out_dir, '{}_{}'.format(i, os.path.basename(bam)))
                isoseq_bams.append((bam, new_bam))
            yield self.clone(IsoSeqBamTrack, trackdb_path=os.path.join(out_dir, 'isoseq_bams.txt'),
                             isoseq_bams=tuple(isoseq_bams))

        if self.genome in pipeline_args.rnaseq_genomes:
            yield self.clone(SpliceTrack, track_path=os.path.join(out_dir, 'splices.bb'),
                             trackdb_path=os.path.join(out_dir, 'splices.txt'))
            # expression is disabled until I fix wiggletools (bamCoverage is needed)
            #if self.genome not in pipeline_args.intron_only_genomes:
            #    yield self.clone(ExpressionTracks, max_track_path=os.path.join(out_dir, 'max_expression.bw'),
            #                     median_track_path=os.path.join(out_dir, 'median_expression.bw'),
            #                     trackdb_path=os.path.join(out_dir, 'expression.txt'))


class CreateTrackDbs(RebuildableTask):
    """Create the final trackDb entries"""
    genome = luigi.Parameter()

    def requires(self):
        return self.clone(CreateTracksDriverTask)

    def output(self):
        pipeline_args = self.get_pipeline_args()
        directory_args = CreateDirectoryStructure.get_args(pipeline_args)
        return luigi.LocalTarget(directory_args.trackdbs[self.genome])

    def run(self):
        pipeline_args = self.get_pipeline_args()
        directory_args = CreateDirectoryStructure.get_args(pipeline_args)
        out_dir = os.path.join(directory_args.out_dir, self.genome)
        org_str = construct_org_str(directory_args.genomes)
        with self.output().open('w') as outf:
            for f in os.listdir(out_dir):
                if f.endswith('.txt'):
                    outf.write('include {}\n'.format(f))
            outf.write('\n\n')

            outf.write(snake_composite.format(org_str=org_str))
            for genome in directory_args.genomes:
                # by default, only the reference genome is visible unless we are on the reference, then all are
                if self.genome == pipeline_args.ref_genome:
                    if genome == pipeline_args.ref_genome:
                        visibility = 'hide'
                    else:
                        visibility = 'full'
                else:
                    visibility = 'hide' if genome != pipeline_args.ref_genome else 'full'
                hal_path = '../{}'.format(os.path.basename(pipeline_args.hal))
                outf.write(snake_template.format(genome=genome, hal_path=hal_path, visibility=visibility))


class DenovoTrack(TrackTask):
    """Constructs a denovo track"""
    mode = luigi.Parameter()

    def run(self):
        def find_rgb(s):
            if s.AssignedGeneId is None and s.AlternativeGeneIds is None:
                return '175,87,207'  # both null -> purple (denovo)
            elif s.AssignedGeneId is None and s.AlternativeGeneIds is not None:
                return '87,207,175'  # no assigned -> teal (possible_paralog)
            return '0'

        def find_alternative_gene_names(s, annotation_info):
            if s.AlternativeGeneIds is None:
                return 'N/A'
            r = {tools.misc.slice_df(annotation_info, gene).iloc[0].GeneName for
                 gene in s.AlternativeGeneIds.split(',')}
            return ','.join(r)

        pipeline_args = self.get_pipeline_args()
        track, trackdb = self.output()
        chrom_sizes = GenomeFiles.get_args(pipeline_args, self.genome).sizes
        # load database information
        db_path = pipeline_args.dbs[self.genome]
        alt_names = load_alt_names(db_path, [self.mode])
        denovo_hgm_df = load_hgm_vectors(db_path, self.mode).drop(['GeneId', 'TranscriptId'], axis=1)
        denovo_df = pd.merge(denovo_hgm_df, alt_names, on='AlignmentId').set_index('AlignmentId')
        annotation_info = tools.sqlInterface.load_annotation(pipeline_args.dbs[pipeline_args.ref_genome])
        annotation_info = annotation_info.set_index('GeneId')

        if self.mode == 'augCGP':
            augustus_gp = AugustusCgp.get_args(pipeline_args).augustus_cgp_gp[self.genome]
        else:
            augustus_gp = AugustusPb.get_args(pipeline_args, self.genome).augustus_pb_gp

        tmp = luigi.LocalTarget(is_tmp=True)
        as_file = luigi.LocalTarget(is_tmp=True)

        with as_file.open('w') as outf:
            outf.write(denovo_as)

        with tmp.open('w') as outf:
            for tx in tools.transcripts.gene_pred_iterator(augustus_gp):
                s = denovo_df.ix[tx.name]
                alternative_gene_ids = 'N/A' if s.AlternativeGeneIds is None else s.AlternativeGeneIds
                intron_rna = ','.join(map(str, s.IntronRnaSupport))
                exon_rna = ','.join(map(str, s.ExonRnaSupport))
                intron_annot = ','.join(map(str, s.IntronAnnotSupport))
                exon_annot = ','.join(map(str, s.ExonAnnotSupport))
                if s.AssignedGeneId is None:
                    assigned_gene_id = gene_name = gene_type = alternative_gene_names = 'N/A'
                else:
                    a = tools.misc.slice_df(annotation_info, s.AssignedGeneId).iloc[0]
                    gene_name = a.GeneName
                    gene_type = a.GeneBiotype
                    assigned_gene_id = s.AssignedGeneId
                    alternative_gene_names = find_alternative_gene_names(s, annotation_info)
                block_starts, block_sizes, exon_frames = tools.transcripts.create_bed_info_gp(tx)
                row = [tx.chromosome, tx.start, tx.stop, tx.name, tx.score, tx.strand, tx.thick_start,
                       tx.thick_stop, find_rgb(s), tx.block_count, block_sizes, block_starts,
                       gene_name, tx.cds_start_stat, tx.cds_end_stat, exon_frames,
                       gene_type, assigned_gene_id, alternative_gene_ids, alternative_gene_names,
                       exon_annot, exon_rna, intron_annot, intron_rna]
                tools.fileOps.print_row(outf, row)
        tools.procOps.run_proc(['bedSort', tmp.path, tmp.path])

        with tools.fileOps.TemporaryFilePath() as out_path:
            cmd = ['bedToBigBed', '-extraIndex=assignedGeneId,name,name2',
                   '-type=bed12+8', '-tab', '-as={}'.format(as_file.path), tmp.path, chrom_sizes, out_path]
            tools.procOps.run_proc(cmd, stderr='/dev/null')
            tools.fileOps.atomic_install(out_path, track.path)

        label = 'AugustusCGP' if self.mode == 'augCGP' else 'AugustusPB'
        description = 'Comparative Augustus' if self.mode == 'augCGP' else 'PacBio Augustus'
        with trackdb.open('w') as outf:
            outf.write(denovo_template.format(name='augustus_{}_{}'.format(self.mode, self.genome),
                                              short_label=label, long_label=label, description=description,
                                              path=os.path.basename(track.path)))


class BgpTrack(TrackTask):
    """Constructs a standard modified bigGenePred track"""
    genepred_path = luigi.Parameter()
    label = luigi.Parameter()
    visibility = luigi.Parameter(default='pack')

    def run(self):
        def find_rgb(info):
            """blue for coding, green for non-coding"""
            if info.TranscriptBiotype == 'protein_coding':
                return '76,85,212'
            return '85,212,76'

        def convert_case(snake_str):
            components = snake_str.split('_')
            return ''.join(x.title() for x in components)

        pipeline_args = self.get_pipeline_args()
        track, trackdb = self.output()
        chrom_sizes = GenomeFiles.get_args(pipeline_args, self.genome).sizes
        annotation_info = tools.sqlInterface.load_annotation(pipeline_args.dbs[pipeline_args.ref_genome])
        # hacky way to make columns consistent
        if 'transcript_id' in annotation_info.columns:
            annotation_info.columns = [convert_case(c) for c in annotation_info.columns]
        annotation_info = annotation_info.set_index('TranscriptId')

        tmp = luigi.LocalTarget(is_tmp=True)
        as_file = luigi.LocalTarget(is_tmp=True)

        with as_file.open('w') as outf:
            outf.write(modified_bgp_as)

        with tmp.open('w') as outf:
            for tx in tools.transcripts.gene_pred_iterator(self.genepred_path):
                s = annotation_info.ix[tools.nameConversions.strip_alignment_numbers(tx.name)]
                block_starts, block_sizes, exon_frames = tools.transcripts.create_bed_info_gp(tx)
                row = [tx.chromosome, tx.start, tx.stop, s.TranscriptName, tx.score, tx.strand, tx.thick_start,
                       tx.thick_stop, find_rgb(s), tx.block_count, block_sizes, block_starts,
                       s.GeneName, tx.cds_start_stat, tx.cds_end_stat, exon_frames,
                       tx.name, s.GeneId, s.TranscriptBiotype, s.GeneBiotype]
                tools.fileOps.print_row(outf, row)
        tools.procOps.run_proc(['bedSort', tmp.path, tmp.path])

        with tools.fileOps.TemporaryFilePath() as out_path:
            cmd = ['bedToBigBed', '-extraIndex=name,name2,geneId,transcriptId',
                   '-type=bed12+8', '-tab', '-as={}'.format(as_file.path), tmp.path, chrom_sizes, out_path]
            tools.procOps.run_proc(cmd, stderr='/dev/null')
            tools.fileOps.atomic_install(out_path, track.path)


        with trackdb.open('w') as outf:
            sanitized_label = self.label.replace(' ', '_').replace('.', '_')
            outf.write(bgp_template.format(name='{}_{}'.format(sanitized_label, self.genome),
                                           label=self.label, visibility=self.visibility,
                                           path=os.path.basename(track.path)))


class ConsensusTrack(TrackTask):
    """Constructs a modified bigGenePred for consensus gene sets"""
    def run(self):
        def find_rgb(info):
            """red for failed, blue for coding, green for non-coding, purple for denovo"""
            if info.transcript_biotype == 'unknown_likely_coding':
                return '135,76,212'
            elif info.transcript_biotype == 'protein_coding':
                return '76,85,212'
            return '85,212,76'

        pipeline_args = self.get_pipeline_args()
        track, trackdb = self.output()
        chrom_sizes = GenomeFiles.get_args(pipeline_args, self.genome).sizes
        consensus_args = Consensus.get_args(pipeline_args, self.genome)
        consensus_gp_info = pd.read_csv(consensus_args.consensus_gp_info, sep='\t',
                                        header=0, na_filter=False).set_index('transcript_id')

        has_rnaseq = len(pipeline_args.rnaseq_genomes) > 0
        has_pb = 'pacbio_isoform_supported' in consensus_gp_info.columns

        tmp_gp = luigi.LocalTarget(is_tmp=True)
        as_file = luigi.LocalTarget(is_tmp=True)

        with tmp_gp.open('w') as outf:
            for tx in tools.transcripts.gene_pred_iterator(consensus_args.consensus_gp):
                info = consensus_gp_info.ix[tx.name]
                block_starts, block_sizes, exon_frames = tools.transcripts.create_bed_info_gp(tx)
                tx_name = info.source_transcript_name if info.source_transcript_name != 'N/A' else tx.name
                row = [tx.chromosome, tx.start, tx.stop, tx_name, tx.score, tx.strand,
                       tx.thick_start, tx.thick_stop, find_rgb(info), tx.block_count, block_sizes, block_starts,
                       info.source_gene_common_name, tx.cds_start_stat, tx.cds_end_stat, exon_frames,
                       tx.name, info.transcript_biotype, tx.name2, info.gene_biotype, info.source_gene,
                       info.source_transcript, info.alignment_id, info.alternative_source_transcripts,
                       info.paralogy, info.frameshift, info.exon_annotation_support,
                       info.intron_annotation_support, info.transcript_class, info.transcript_modes,
                       info.valid_start, info.valid_stop, info.proper_orf]
                if has_rnaseq:
                    row.extend([info.intron_rna_support, info.exon_rna_support])
                if has_pb:
                    row.append(info.pacbio_isoform_supported)
                tools.fileOps.print_row(outf, row)

        with as_file.open('w') as outf:
            as_str = construct_consensus_gp_as(has_rnaseq, has_pb)
            outf.write(as_str)
        tools.procOps.run_proc(['bedSort', tmp_gp.path, tmp_gp.path])
        with tools.fileOps.TemporaryFilePath() as out_path:
            cmd = ['bedToBigBed', '-extraIndex=name,name2,txId,geneName,sourceGene,sourceTranscript,alignmentId',
                   '-type=bed12+21', '-tab', '-as={}'.format(as_file.path), tmp_gp.path, chrom_sizes, out_path]
            tools.procOps.run_proc(cmd, stderr='/dev/null')
            tools.fileOps.atomic_install(out_path, track.path)

        with trackdb.open('w') as outf:
            outf.write(consensus_template.format(genome=self.genome, path=os.path.basename(track.path)))


class EvaluationTrack(TrackTask):
    """Constructs the consensus evaluation track"""
    tx_modes = luigi.TupleParameter()

    def run(self):
        def load_evals(tx_mode):
            """Loads the error tracks from the database"""
            cds_table = tools.sqlInterface.tables['CDS'][tx_mode]['evaluation']
            mrna_table = tools.sqlInterface.tables['mRNA'][tx_mode]['evaluation']
            cds_df = pd.read_sql_table(cds_table.__tablename__, engine).set_index('AlignmentId')
            mrna_df = pd.read_sql_table(mrna_table.__tablename__, engine).set_index('AlignmentId')
            return {'CDS': cds_df, 'mRNA': mrna_df}

        pipeline_args = self.get_pipeline_args()
        track, trackdb = self.output()
        chrom_sizes = GenomeFiles.get_args(pipeline_args, self.genome).sizes
        engine = tools.sqlInterface.create_engine('sqlite:///' + pipeline_args.dbs[self.genome])
        evals = {tx_mode: load_evals(tx_mode) for tx_mode in self.tx_modes}
        consensus_args = Consensus.get_args(pipeline_args, self.genome)
        consensus_gp_info = pd.read_csv(consensus_args.consensus_gp_info, sep='\t',
                                        header=0, na_filter=False).set_index('transcript_id')
        aln_ids = set(consensus_gp_info.alignment_id)
        rows = []
        for aln_id in aln_ids:
            tx_mode = tools.nameConversions.alignment_type(aln_id)
            if tx_mode not in ['transMap', 'augTM', 'augTMR']:
                continue
            mode = 'CDS'
            df = tools.misc.slice_df(evals[tx_mode][mode], aln_id)
            if len(df) == 0:
                mode = 'mRNA'
                df = tools.misc.slice_df(evals[tx_mode][mode], aln_id)
            for tx_id, s in df.iterrows():
                bed = s.tolist()
                bed[3] = '/'.join([tx_id, bed[3], mode])
                rows.append(bed)

        tmp = luigi.LocalTarget(is_tmp=True)
        with tmp.open('w') as tmp_handle:
            tools.fileOps.print_rows(tmp_handle, rows)
        tools.procOps.run_proc(['bedSort', tmp.path, tmp.path])
        with tools.fileOps.TemporaryFilePath() as out_path:
            cmd = ['bedToBigBed', '-type=bed12', '-tab', tmp.path, chrom_sizes, out_path]
            tools.procOps.run_proc(cmd, stderr='/dev/null')
            tools.fileOps.atomic_install(out_path, track.path)


        with trackdb.open('w') as outf:
            outf.write(error_template.format(genome=self.genome, path=os.path.basename(track.path)))


class TransMapTrack(TrackTask):
    """Constructs the transMap bigPsl"""
    def run(self):
        pipeline_args = self.get_pipeline_args()
        track, trackdb = self.output()
        chrom_sizes = GenomeFiles.get_args(pipeline_args, self.genome).sizes
        fasta = ReferenceFiles.get_args(pipeline_args).transcript_fasta
        tm_args = TransMap.get_args(pipeline_args, self.genome)

        # we have to construct a reference CDS file as well as rename the reference fasta to the new names
        mrna = luigi.LocalTarget(is_tmp=True)
        cds = luigi.LocalTarget(is_tmp=True)
        tmp = luigi.LocalTarget(is_tmp=True)
        as_file = luigi.LocalTarget(is_tmp=True)
        seq_dict = tools.bio.get_sequence_dict(fasta)
        ref_tx_dict = tools.transcripts.get_gene_pred_dict(ReferenceFiles.get_args(pipeline_args).annotation_gp)
        with cds.open('w') as cds_handle, mrna.open('w') as mrna_handle:
            for tx in tools.transcripts.gene_pred_iterator(tm_args.tm_gp):
                ref_tx = ref_tx_dict[tools.nameConversions.strip_alignment_numbers(tx.name)]
                tools.bio.write_fasta(mrna_handle, tx.name, str(seq_dict[ref_tx.name]))
                if ref_tx.cds_size == 0:  # non-coding txs have no cds interval
                    start = stop = 0
                else:
                    start = ref_tx.cds_coordinate_to_mrna(0) + ref_tx.offset
                    stop = start + ref_tx.cds_size - ((ref_tx.cds_size - ref_tx.offset) % 3)
                cds_handle.write('{}\t{}..{}\n'.format(tx.name, start + 1, stop))

        with tmp.open('w') as outf:
            cmd = [['pslToBigPsl', '-cds={}'.format(cds.path), '-fa={}'.format(mrna.path), tm_args.tm_psl, 'stdout'],
                   ['bedSort', '/dev/stdin', '/dev/stdout']]
            tools.procOps.run_proc(cmd, stdout=outf, stderr='/dev/null')

        with as_file.open('w') as outf:
            outf.write(bigpsl)

        with tools.fileOps.TemporaryFilePath() as out_path:
            cmd = ['bedToBigBed', '-type=bed12+13', '-tab', '-extraIndex=name',
                   '-as={}'.format(as_file.path), tmp.path, chrom_sizes, out_path]
            tools.procOps.run_proc(cmd, stderr='/dev/null')
            tools.fileOps.atomic_install(out_path, track.path)

        with trackdb.open('w') as outf:
            outf.write(bigpsl_template.format(name='transmap_{}'.format(self.genome), short_label='transMap',
                                              long_label='transMap', path=os.path.basename(track.path),
                                              visibility='pack'))


class AugustusTrack(TrackTask):
    """Constructs a combined TM(R) track"""
    def run(self):
        pipeline_args = self.get_pipeline_args()
        track, trackdb = self.output()
        chrom_sizes = GenomeFiles.get_args(pipeline_args, self.genome).sizes
        annotation_info = tools.sqlInterface.load_annotation(pipeline_args.dbs[pipeline_args.ref_genome])
        annotation_info = annotation_info.set_index('TranscriptId')
        aug_args = Augustus.get_args(pipeline_args, self.genome)
        tm_gp = aug_args.augustus_tm_gp
        if self.genome in pipeline_args.rnaseq_genomes:
            tmr_gp = aug_args.augustus_tmr_gp
        else:
            tmr_gp = None

        with tools.fileOps.TemporaryFilePath() as tmp, tools.fileOps.TemporaryFilePath() as as_file:
            with open(as_file, 'w') as outf:
                outf.write(modified_bgp_as)
            with open(tmp, 'w') as outf:
                for gp, color in zip(*[[tm_gp, tmr_gp], ['38,112,75', '112,38,75']]):
                    if gp is None:
                        continue
                    gp = tools.transcripts.gene_pred_iterator(gp)
                    for tx in gp:
                        s = annotation_info.ix[tools.nameConversions.strip_alignment_numbers(tx.name)]
                        block_starts, block_sizes, exon_frames = tools.transcripts.create_bed_info_gp(tx)
                        row = [tx.chromosome, tx.start, tx.stop, s.TranscriptName, tx.score, tx.strand, tx.thick_start,
                               tx.thick_stop, color, tx.block_count, block_sizes, block_starts,
                               s.GeneName, tx.cds_start_stat, tx.cds_end_stat, exon_frames,
                               tx.name, s.GeneId, s.TranscriptBiotype, s.GeneBiotype]
                        tools.fileOps.print_row(outf, row)
            tools.procOps.run_proc(['bedSort', tmp, tmp])

            with tools.fileOps.TemporaryFilePath() as out_path:
                cmd = ['bedToBigBed', '-extraIndex=name,name2,geneId,transcriptId',
                       '-type=bed12+8', '-tab', '-as={}'.format(as_file), tmp, chrom_sizes, out_path]
                tools.procOps.run_proc(cmd, stderr='/dev/null')
                tools.fileOps.atomic_install(out_path, track.path)

        with trackdb.open('w') as outf:
            outf.write(bgp_template.format(name='augustus_{}'.format(self.genome), label='AugustusTM(R)',
                                           path=os.path.basename(track.path), visibility='hide'))


class IsoSeqBamTrack(RebuildableTask):
    """Symlinks over IsoSeq bams"""
    genome = luigi.Parameter()
    trackdb_path = luigi.Parameter()
    isoseq_bams = luigi.TupleParameter()

    def output(self):
        r = [luigi.LocalTarget(new_bam) for bam, new_bam in self.isoseq_bams]
        r.extend([luigi.LocalTarget(x.path + '.bai') for x in r])
        return r, luigi.LocalTarget(self.trackdb_path)

    def requires(self):
        return self.clone(CreateDirectoryStructure), self.clone(Consensus)

    def run(self):
        bams, trackdb = self.output()
        with trackdb.open('w') as outf:
            outf.write(bam_composite_template.format(genome=self.genome))
            for bam, new_bam in self.isoseq_bams:
                shutil.copy(bam, new_bam)
                shutil.copy(bam + '.bai', new_bam + '.bai')
                name = os.path.splitext(os.path.basename(bam))[0].split('_', 1)[0]
                outf.write(bam_template.format(bam=os.path.basename(new_bam), name=name, genome=self.genome))


class SpliceTrack(TrackTask):
    """Constructs the splice junction track"""
    def run(self):
        def parse_entry(entry):
            """Converts a GFF entry to BED12"""
            start = int(entry[3]) - 1
            stop = int(entry[4])
            block_starts = '0,{}'.format(stop - start - 2)
            mult = int(tools.misc.parse_gff_attr_line(entry[-1])['mult'])
            return [entry[0], start, stop, 'SpliceJunction', mult, '.', start, stop, '204,124,45',
                    '2', '2,2', block_starts]

        pipeline_args = self.get_pipeline_args()
        track, trackdb = self.output()
        chrom_sizes = GenomeFiles.get_args(pipeline_args, self.genome).sizes
        hints_gff = BuildDb.get_args(pipeline_args, self.genome).hints_path

        entries = []
        for line in open(hints_gff):
            if '\tintron\t' in line and 'src=E' in line and 'mult' in line:
                parsed = parse_entry(line.split('\t'))
                if parsed[4] > 2:
                    entries.append(parsed)

        mults = [x[4] for x in entries]
        tot = sum(mults)
        # calculate junctions per thousand
        jpt = [1.0 * x * 10 ** 4 / tot for x in mults]
        for e, s in zip(*[entries, jpt]):
            e[4] = min(int(round(s)), 1000)

        # load to file
        tmp = luigi.LocalTarget(is_tmp=True)
        with tmp.open('w') as tmp_handle:
            tools.fileOps.print_rows(tmp_handle, entries)

        tools.procOps.run_proc(['bedSort', tmp.path, tmp.path])

        with tools.fileOps.TemporaryFilePath() as out_path:
            cmd = ['bedToBigBed', '-tab', tmp.path, chrom_sizes, out_path]
            tools.procOps.run_proc(cmd, stderr='/dev/null')
            tools.fileOps.atomic_install(out_path, track.path)

        with trackdb.open('w') as outf:
            outf.write(splice_template.format(genome=self.genome, path=os.path.basename(track.path)))


class ExpressionTracks(RebuildableTask):
    """Constructs the maximum and median expression tracks"""
    genome = luigi.Parameter()
    trackdb_path = luigi.Parameter()
    max_track_path = luigi.Parameter()
    median_track_path = luigi.Parameter()

    def output(self):
        return [luigi.LocalTarget(self.max_track_path), luigi.LocalTarget(self.median_track_path)], \
               luigi.LocalTarget(self.trackdb_path)

    def requires(self):
        return self.clone(CreateDirectoryStructure), self.clone(Consensus)

    def run(self):
        pipeline_args = self.get_pipeline_args()
        bams = list(pipeline_args.cfg['BAM'][self.genome])
        (max_track, median_track), trackdb = self.output()
        chrom_sizes = GenomeFiles.get_args(pipeline_args, self.genome).sizes

        with median_track.open('w') as outf:
            cmd = [['wiggletools', 'median'] + bams,
                   ['wigToBigWig', '-clip', '/dev/stdin', chrom_sizes, '/dev/stdout']]
            tools.procOps.run_proc(cmd, stdout=outf, stderr='/dev/null')

        with max_track.open('w') as outf:
            cmd = [['wiggletools', 'max'] + bams,
                   ['wigToBigWig', '-clip', '/dev/stdin', chrom_sizes, '/dev/stdout']]
            tools.procOps.run_proc(cmd, stdout=outf, stderr='/dev/null')

        with trackdb.open('w') as outf:
            outf.write(wiggle_template.format(genome=self.genome, mode='Median',
                                              path=os.path.basename(median_track.path), color='151,189,68'))
            outf.write(wiggle_template.format(genome=self.genome, mode='Maximum',
                                              path=os.path.basename(max_track.path), color='106,68,189'))


###
# assembly hub functions, including generating autoSql files
###


def find_default_pos(chrom_sizes, window_size=200000):
    """
    Returns a window_size window over the beginning of the largest chromosome
    :param chrom_sizes: chrom sizes file
    :param window_size: window size to extend from
    :return: string
    """
    sizes = [x.split() for x in open(chrom_sizes)]
    sorted_sizes = sorted(sizes, key=lambda (chrom, size): -int(size))
    return '{}:{}-{}'.format(sorted_sizes[0][0], 1, window_size)


def construct_org_str(genomes):
    """Constructs the organism string for the hal snakes. format is genome=genome space separated"""
    return ' '.join(['{0}={0}'.format(genome) for genome in genomes])


def construct_consensus_gp_as(has_rna, has_pb):
    """Dynamically generate an autosql file for consensus"""
    consensus_gp_as = '''table bigCat
"bigCat gene models"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name"
    uint score;         "Score (0-1000)"
    char[1] strand;     "+ or - for strand"
    uint thickStart;    "Start of where display should be thick (start codon)"
    uint thickEnd;      "End of where display should be thick (stop codon)"
    uint reserved;       "RGB value (use R,G,B string in input file)"
    int blockCount;     "Number of blocks"
    int[blockCount] blockSizes; "Comma separated list of block sizes"
    int[blockCount] chromStarts; "Start positions relative to chromStart"
    string name2;       "Gene name"
    string cdsStartStat; "Status of CDS start annotation"
    string cdsEndStat;   "Status of CDS end annotation"
    int[blockCount] exonFrames; "Exon frame {0,1,2}, or -1 if no frame for exon"
    string txId; "Transcript ID"
    string type;        "Transcript type"
    string geneName;    "Gene ID"
    string geneType;    "Gene type"
    string sourceGene;    "Source gene ID"
    string sourceTranscript;    "Source transcript ID"
    string alignmentId;  "Alignment ID"
    lstring alternativeSourceTranscripts;    "Alternative source transcripts"
    lstring Paralogy;    "Paralogous alignment IDs"
    string frameshift;  "Frameshifted relative to source?"
    lstring exonAnnotationSupport;   "Exon support in reference annotation"
    lstring intronAnnotationSupport;   "Intron support in reference annotation"
    string transcriptClass;    "Transcript class"
    string transcriptModes;    "Transcript mode(s)"
    string validStart;         "Valid start codon"
    string validStop;          "Valid stop codon"
    string properOrf;           "Proper multiple of 3 ORF"
'''
    if has_rna:
        consensus_gp_as += '    lstring intronRnaSupport;   "RNA intron support"\n'
        consensus_gp_as += '    lstring exonRnaSupport;  "RNA exon support"\n'
    if has_pb:
        consensus_gp_as += '    string pbIsoformSupported;   "Is this transcript supported by IsoSeq?"'
    consensus_gp_as += '\n)\n'
    return consensus_gp_as


modified_bgp_as = '''table bigGenePred
"bigGenePred gene models"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name"
    uint score;         "Score (0-1000)"
    char[1] strand;     "+ or - for strand"
    uint thickStart;    "Start of where display should be thick (start codon)"
    uint thickEnd;      "End of where display should be thick (stop codon)"
    uint reserved;       "RGB value (use R,G,B string in input file)"
    int blockCount;     "Number of blocks"
    int[blockCount] blockSizes; "Comma separated list of block sizes"
    int[blockCount] chromStarts; "Start positions relative to chromStart"
    string name2;       "Gene name"
    string cdsStartStat; "Status of CDS start annotation (none, unknown, incomplete, or complete)"
    string cdsEndStat;   "Status of CDS end annotation (none, unknown, incomplete, or complete)"
    int[blockCount] exonFrames; "Exon frame {0,1,2}, or -1 if no frame for exon"
    string transcriptId;  "Transcript ID"
    string geneId;    "Gene ID"
    string type;        "Transcript type"
    string geneType;    "Gene type"
    )

'''


denovo_as = '''table denovo
"denovo gene models"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name"
    uint score;         "Score (0-1000)"
    char[1] strand;     "+ or - for strand"
    uint thickStart;    "Start of where display should be thick (start codon)"
    uint thickEnd;      "End of where display should be thick (stop codon)"
    uint reserved;       "RGB value (use R,G,B string in input file)"
    int blockCount;     "Number of blocks"
    int[blockCount] blockSizes; "Comma separated list of block sizes"
    int[blockCount] chromStarts; "Start positions relative to chromStart"
    string name2;       "Assigned gene name"
    string cdsStartStat; "Status of CDS start annotation (none, unknown, incomplete, or complete)"
    string cdsEndStat;   "Status of CDS end annotation (none, unknown, incomplete, or complete)"
    int[blockCount] exonFrames; "Exon frame {0,1,2}, or -1 if no frame for exon"
    string geneType;    "Assigned gene type"
    string assignedGeneId; "Assigned source gene ID"
    lstring alternativeGeneIds; "Alternative source gene IDs"
    lstring alternativeGeneNames; "Alternative source gene names"
    lstring exonAnnotationSupport;   "Exon support in reference annotation"
    lstring exonRnaSupport;  "RNA exon support"
    lstring intronAnnotationSupport;   "Intron support in reference annotation"
    lstring intronRnaSupport;   "RNA intron support"
    )
'''


bigpsl = '''table bigPsl
"bigPsl pairwise alignment"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name"
    uint score;         "Score (0-1000)"
    char[1] strand;     "+ or - indicates whether the query aligns to the + or - strand on the reference"
    uint thickStart;    "Start of where display should be thick (start codon)"
    uint thickEnd;      "End of where display should be thick (stop codon)"
    uint reserved;       "RGB value (use R,G,B string in input file)"
    int blockCount;     "Number of blocks"
    int[blockCount] blockSizes; "Comma separated list of block sizes"
    int[blockCount] chromStarts; "Start positions relative to chromStart"
    uint    oChromStart;"Start position in other chromosome"
    uint    oChromEnd;  "End position in other chromosome"
    char[1] oStrand;    "+ or -, - means that psl was reversed into BED-compatible coordinates"
    uint    oChromSize; "Size of other chromosome."
    int[blockCount] oChromStarts; "Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand"
    lstring  oSequence;  "Sequence on other chrom (or edit list, or empty)"
    string   oCDS;       "CDS in NCBI format"
    uint    chromSize;"Size of target chromosome"
    uint match;        "Number of bases matched."
    uint misMatch; " Number of bases that don't match "
    uint repMatch; " Number of bases that match but are part of repeats "
    uint nCount;   " Number of 'N' bases "
    uint seqType;    "0=empty, 1=nucleotide, 2=amino_acid"
    )

'''


###
# Templates for trackDb entries
###


hub_str = '''hub {hal}
shortLabel {hal}
longLabel {hal}
genomesFile genomes.txt
email {email}

'''

genome_str = '''genome {genome}
twoBitPath {genome}/{genome}.2bit
trackDb {genome}/trackDb.txt
organism {genome}
description {genome}
scientificName {genome}
defaultPos {default_pos}
groups groups.txt

'''


groups_str = '''name cat_tracks
label Comparative Annotation Toolkit
priority 1
defaultIsClosed 0

name snake
label Alignment Snakes
priority 2
defaultIsClosed 0

name expression
label Expression
priority 3
defaultIsClosed 0

'''


snake_composite = '''track hubCentral
compositeTrack on
shortLabel Cactus
longLabel Cactus Alignment Tracks
group cat_tracks
subGroup1 view Track_Type Snake=Alignments
subGroup2 orgs Organisms {org_str}
dragAndDrop subTracks
dimensions dimensionX=view dimensionY=orgs
noInherit on
priority 0
centerLabelsDense on
visibility full
type bigBed 3

    track hubCentralAlignments
    shortLabel Alignments
    view Alignments
    visibility full
    subTrack hubCentral

'''

snake_template = '''        track snake{genome}
        longLabel {genome}
        shortLabel {genome}
        otherSpecies {genome}
        visibility {visibility}
        parent hubCentralAlignments off
        priority 3
        bigDataUrl {hal_path}
        type halSnake
        group snake
        subGroups view=Snake orgs={genome}

'''


consensus_template = '''track consensus_{genome}
shortLabel CAT Annotation
longLabel CAT Annotation
description CAT Annotation
type bigGenePred
group cat_tracks
itemRgb on
priority 1
visibility pack
searchIndex name,name2,txId,geneName,sourceGene,sourceTranscript,alignmentId
bigDataUrl {path}
labelFields name,name2,txId,geneName,sourceGene,sourceTranscript,alignmentId
defaultLabelFields name
labelSeperator " "

'''


bgp_template = '''track {name}
shortLabel {label}
longLabel {label}
description {label}
type bigGenePred
group cat_tracks
itemRgb on
priority 3
visibility {visibility}
searchIndex name,name2,geneId,transcriptId
bigDataUrl {path}
labelFields name,name2,geneId,transcriptId
defaultLabelFields name
labelSeperator " "

'''


bigpsl_template = '''track {name}
shortLabel {short_label}
longLabel {long_label}
bigDataUrl {path}
type bigPsl
group cat_tracks
itemRgb on
priority 2
visibility {visibility}
baseColorUseSequence lfExtra
baseColorDefault genomicCodons
baseColorUseCds given
indelDoubleInsert on
indelQueryInsert on
showDiffBasesAllScales .
showDiffBasesMaxZoom 10000.0
#showCdsAllScales .
#showCdsMaxZoom 10000.0
searchIndex name

'''


denovo_template = '''track {name}
shortLabel {short_label}
longLabel {long_label}
description {description}
bigDataUrl {path}
type bigGenePred
group cat_tracks
priority 4
itemRgb on
searchIndex assignedGeneId,name,name2
labelFields assignedGeneId,name,name2
defaultLabelFields name
labelSeperator " "


'''


bam_composite_template = '''track bams_{genome}
group expression
compositeTrack on
shortLabel IsoSeq BAMs
longLabel IsoSeq BAMs
dragAndDrop subTracks
visibility hide
type bam
indelDoubleInsert on
indelQueryInsert on
showNames off
bamColorMode gray
bamGrayMode aliQual
pairEndsByName on

'''

bam_template = '''    track {bam}_{genome}
    parent bams_{genome}
    bigDataUrl {bam}
    shortLabel {name}
    longLabel {name}
    type bam
    priority 10

'''


wiggle_template = '''track {mode}_{genome}
shortLabel {mode} expression
longLabel {mode} expression
type bigWig
group expression
bigDataUrl {path}
color {color}
visibility hide
priority 11
spectrum on

'''


splice_template = '''track splices_{genome}
type bigBed 12
group expression
shortLabel RNA-seq splices
longLabel RNA-seq Splice Junctions
bigDataUrl {path}
visibility hide
color 45,125,204
priority 12
spectrum on
minGrayLevel 4

'''


error_template = '''track error_{genome}
type bigBed 6
group cat_tracks
shortLabel Consensus indels
longLabel Consensus indels
bigDataUrl {path}
visibility hide
priority 5

'''
