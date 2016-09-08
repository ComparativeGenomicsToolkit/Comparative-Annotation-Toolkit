"""
Addons for Luigi. Includes decorators multiple inheritance and requirements as well as abstract classes extending
both the Task and Target paradigms.
"""
import os
import hal
import argparse
import luigi
import luigi.util
import luigi.contrib.sqla
import procOps


class PipelineTask(luigi.Task):
    """
    Base class for all tasks in this pipeline. Provides useful methods for handling parameters being passed between
    modules.
    """

    def get_pipeline_args(self):
        """returns a args of all of the arguments to the pipeline. Resolves the target genomes variable"""
        args = argparse.Namespace()
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
        if self.augustus_cgp_cfg is not None:
            args.augustus_cgp_cfg = os.path.abspath(self.augustus_cgp_cfg)
        else:
            args.augustus_cgp_cfg = None
        if self.augustus_cgp_param is not None:
            args.augustus_cgp_param = os.path.abspath(self.augustus_cgp_param)
        else:
            args.augustus_cgp_param = None
        args.hal_genomes = hal.extract_genomes(self.hal)
        if self.target_genomes is None:
            target_genomes = tuple(set(args.hal_genomes) - {self.ref_genome})
        else:
            target_genomes = tuple([x for x in self.target_genomes])
        args.target_genomes = target_genomes
        args.modes = self.get_modes(args)
        args.dbs = self.get_databases(args)
        return args

    def get_modes(self, pipeline_args):
        """Convenience function that reports the modes we are operating in as a list"""
        modes = ['transMap']
        if pipeline_args.augustus is True:
            modes.append('augTM')
        if pipeline_args.augustus_tmr is True:
            modes.append('augTMR')
        if pipeline_args.augustus_cgp is True:
            modes.append('augCGP')
        return modes

    def get_databases(self, pipeline_args):
        """wrapper for get_database() that provides all of the databases"""
        dbs = {genome: self.get_database(pipeline_args, genome) for genome in pipeline_args.hal_genomes}
        return dbs

    def get_database(self, pipeline_args, genome):
        """database paths must be resolved here to handle multiple programs accessing them"""
        base_out_dir = os.path.join(pipeline_args.out_dir, 'databases')
        return os.path.join(base_out_dir, '{}.db'.format(genome))

    def get_module_args(self, module, **args):
        """
        convenience wrapper that takes a parent module and propagates any required arguments to generate the full
        argument set.
        """
        pipeline_args = self.get_pipeline_args()
        return module.get_args(pipeline_args, **args)


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
            procOps.run_proc(cmd, stdout=outf)


class multiple_inherits(object):
    """
    Task inheritance.

    Usage:

    .. code-block:: python

        class TaskPathA(luigi.Task):
            a = luigi.IntParameter()
            # ...

        class TaskPathB(luigi.Task):
            b = luigi.IntParameter()

        @multiple_inherits(TaskPathA, TaskPathB):
        class MyTask(luigi.Task):
            def requires(self):
               return self.clone_parent()

            def run(self):
               print self.a # this will be defined
               print self.b # this will also be defined
               # ...
    """

    def __init__(self, *tasks_to_inherit):
        super(multiple_inherits, self).__init__()
        self.tasks_to_inherit = tasks_to_inherit

    def __call__(self, task_that_inherits):
        tasks_to_inherit = self.tasks_to_inherit
        for task_to_inherit in tasks_to_inherit:
            for param_name, param_obj in task_to_inherit.get_params():
                if not hasattr(task_that_inherits, param_name):
                    setattr(task_that_inherits, param_name, param_obj)

        # Modify task_that_inherits by subclassing it and adding methods
        @luigi.util.task_wraps(task_that_inherits)
        class Wrapped(task_that_inherits):
            def clone_parent(self, **args):
                task = self.clone(cls=tasks_to_inherit[0])
                for additional_task in tasks_to_inherit[1:]:
                    task = task.clone(cls=additional_task, **args)
                return task
        return Wrapped


class multiple_requires(object):
    """
    Same as @multiple_inherits, but also auto-defines the requires method.
    """

    def __init__(self, *tasks_to_require):
        super(multiple_requires, self).__init__()
        self.inherit_decorator = multiple_inherits(*tasks_to_require)
        self.tasks_to_require = tasks_to_require

    def __call__(self, task_that_requires):
        task_that_requires = self.inherit_decorator(task_that_requires)
        tasks_to_require = self.tasks_to_require

        # Modify task_that_requires by subclassing it and adding methods
        @luigi.util.task_wraps(task_that_requires)
        class Wrapped(task_that_requires):
            def requires(self):
                return (self.clone(x) for x in tasks_to_require)
        return Wrapped
