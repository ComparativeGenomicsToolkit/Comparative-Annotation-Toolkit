"""
Base tasks used by the pipeline. 
"""
import argparse
import json
import os
import shutil
import tempfile

import luigi
from frozendict import frozendict
from toil.job import Job

import tools.fileOps
import tools.procOps


class HashableNamespace(argparse.Namespace):
    """
    Adds a __hash__ function to argparse's Namespace. Follows best practices for implementation of __hash__.
    """
    def __hash__(self):
        def xor(x, y):
            return x ^ hash(y)
        val_iter = self.__dict__.itervalues()
        first = hash(val_iter.next())
        return reduce(xor, val_iter, first) ^ hash(tuple(self.__dict__.values()))


class PipelineTask(luigi.Task):
    """
    Base class for all tasks in this pipeline. Has Parameters for all input parameters that will be inherited
    by all downstream tools.

    Provides useful methods for handling parameters being passed between modules.
    """
    hal = luigi.Parameter()
    ref_genome = luigi.Parameter()
    annotation = luigi.Parameter()
    out_dir = luigi.Parameter(default='./cat_output')
    work_dir = luigi.Parameter(default=os.path.join(tempfile.gettempdir(), __name__))
    target_genomes = luigi.TupleParameter(default=None)
    # Debugging option - use this to bypass the dependency graph for specific submodules
    no_evaluate_dependency = luigi.BoolParameter(default=False)
    # AugustusTM(R) parameters
    augustus = luigi.BoolParameter(default=False)
    augustus_species = luigi.Parameter(default='human', significant=False)
    augustus_hints_db = luigi.Parameter(default=None)
    tm_cfg = luigi.Parameter(default='augustus_cfgs/extrinsic.ETM1.cfg', significant=False)
    tmr_cfg = luigi.Parameter(default='augustus_cfgs/extrinsic.ETM2.cfg', significant=False)
    # AugustusCGP parameters
    augustus_cgp = luigi.BoolParameter(default=False)
    cgp_param = luigi.Parameter(default='augustus_cfgs/log_reg_parameters_default.cfg', significant=False)
    augustus_cgp_cfg_template = luigi.Parameter(default='augustus_cfgs/cgp_extrinsic_template.cfg', significant=False)
    maf_chunksize = luigi.IntParameter(default=2500000, significant=False)
    maf_overlap = luigi.IntParameter(default=500000, significant=False)
    # consensus options
    resolve_split_genes = luigi.BoolParameter(default=False)
    # Toil options
    batchSystem = luigi.Parameter(default='singleMachine', significant=False)
    maxCores = luigi.IntParameter(default=32, significant=False)
    logLevel = luigi.Parameter(default='WARNING', significant=False)  # this is passed to toil
    cleanWorkDir = luigi.Parameter(default='onSuccess', significant=False)  # debugging option
    parasolCommand = luigi.Parameter(default=None, significant=False)
    defaultMemory = luigi.IntParameter(default=8 * 1024 ** 3, significant=False)
    workDir = luigi.Parameter(default=None, significant=False)
    disableCaching = luigi.BoolParameter(default=False, significant=False)

    def __repr__(self):
        """override the repr to make logging cleaner"""
        # we are in a genome-specific task, so say so
        if hasattr(self, 'genome'):
            return 'Task: {} for {}'.format(self.__class__.__name__, self.genome)
        elif hasattr(self, 'mode'):
            return 'Task: {} for {}'.format(self.__class__.__name__, self.mode)
        else:
            return 'Task: {}'.format(self.__class__.__name__)

    def get_pipeline_args(self):
        """returns a namespace of all of the arguments to the pipeline. Resolves the target genomes variable"""
        args = HashableNamespace()
        args.hal = os.path.abspath(self.hal)
        args.ref_genome = self.ref_genome
        args.annotation = os.path.abspath(self.annotation)
        args.out_dir = os.path.abspath(self.out_dir)
        args.work_dir = os.path.abspath(self.work_dir)
        args.augustus = self.augustus
        args.augustus_cgp = self.augustus_cgp
        args.augustus_species = self.augustus_species
        if self.augustus_hints_db is not None:
            args.augustus_hints_db = os.path.abspath(self.augustus_hints_db)
            args.augustus_tmr = True if self.augustus else False
        else:
            args.augustus_hints_db = None
        args.tm_cfg = os.path.abspath(self.tm_cfg)
        args.tmr_cfg = os.path.abspath(self.tmr_cfg)
        args.augustus_cgp = self.augustus_cgp
        args.maf_chunksize = self.maf_chunksize
        args.maf_overlap = self.maf_overlap
        args.resolve_split_genes = self.resolve_split_genes
        args.augustus_cgp_cfg_template = os.path.abspath(self.augustus_cgp_cfg_template)
        if self.cgp_param is not None:
            args.cgp_param = os.path.abspath(self.cgp_param)
        else:
            args.cgp_param = None
        args.hal_genomes = tuple(tools.hal.extract_genomes(self.hal))
        if self.target_genomes is None:
            target_genomes = tuple(set(args.hal_genomes) - {self.ref_genome})
        else:
            target_genomes = tuple([x for x in self.target_genomes])
        args.target_genomes = target_genomes
        args.modes = PipelineTask.get_modes(args)
        args.dbs = PipelineTask.get_databases(args)
        args.max_cores = self.maxCores  # used for HGM
        return args

    def get_module_args(self, module, **args):
        """
        convenience wrapper that takes a parent module and propagates any required arguments to generate the full
        argument set.
        """
        pipeline_args = self.get_pipeline_args()
        return module.get_args(pipeline_args, **args)

    @staticmethod
    def get_modes(pipeline_args):
        """Convenience function that reports the modes we are operating in as a list"""
        modes = ['transMap']
        if pipeline_args.augustus is True:
            modes.append('augTM')
        if pipeline_args.augustus_tmr is True:
            modes.append('augTMR')
        if pipeline_args.augustus_cgp is True:
            modes.append('augCGP')
        return tuple(modes)

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
        # luigi localTargets guarantee atomicity if used as a handle
        with self.output().open('w') as outf:
            tools.procOps.run_proc(cmd, stdout=outf)


class ToilTask(PipelineTask):
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
        if toil_args.workDir is not None:
            if toil_args.batchSystem == 'parasol' and toil_args.disableCaching is False:
                raise RuntimeError('Running parasol without disabled caching is a very bad idea.')
            tools.fileOps.ensure_dir(toil_args.workDir)
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
        base_repr = super(ToilTask, self).__repr__()
        return 'Toil' + base_repr + ' using batchSystem {}'.format(self.batchSystem)


class HintsDbTask(luigi.Task):
    # path to config file
    config = luigi.Parameter()
    hal = luigi.Parameter()
    augustus_hints_db = luigi.Parameter(default='augustus_hints.db')
    work_dir = luigi.Parameter(default=os.path.join(tempfile.gettempdir(), __name__))
    use_wiggle_hints = luigi.BoolParameter(default=True)
    # Toil options
    batchSystem = luigi.Parameter(default='singleMachine', significant=False)
    maxCores = luigi.IntParameter(default=16, significant=False)
    logLevel = luigi.Parameter(default='WARNING', significant=False)  # this is passed to toil
    cleanWorkDir = luigi.Parameter(default='onSuccess', significant=False)  # debugging option
    workDir = luigi.Parameter(default=None, significant=False)
    disableCaching = luigi.BoolParameter(default=False, significant=False)
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
