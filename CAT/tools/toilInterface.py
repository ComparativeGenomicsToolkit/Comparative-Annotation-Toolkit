"""
Implements a hashable and iterable namespace.
Functions for validation using the ConfigObj library.
"""
import os
import tempfile
import luigi
import tools.fileOps
from toil.job import Job


class ToilTask(luigi.Task):
    """
    Task for launching toil pipelines from within luigi.
    """
    workDir = luigi.Parameter(default=tempfile.gettempdir())
    batchSystem = luigi.Parameter(default='singleMachine')
    maxCores = luigi.IntParameter(default=16)
    logLevel = luigi.Parameter(default='WARNING')

    def prepare_toil_options(self, job_sub_path):
        """
        Prepares a Namespace object for Toil which has all defaults, overridden as specified
        Will see if the jobStore path exists, and if it does, assume that we need to add --restart
        :param job_sub_path: The sub path that will be joined to make a full path. For example, 'chaining/C57BL_6NJ'
        :return: Namespace
        """
        job_store = os.path.abspath(os.path.join(self.workDir, 'toil', job_sub_path))
        toil_args = get_toil_defaults()
        toil_args.__dict__.update(vars(self))
        toil_args.jobStore = job_store
        if os.path.exists(job_store):
            toil_args.restart = True
        else:
            tools.fileOps.ensure_file_dir(job_store)
        return toil_args


def get_toil_defaults():
    """
    Extracts the default toil options as a dictionary, setting jobStore to None
    :return: dict
    """
    parser = Job.Runner.getDefaultArgumentParser()
    namespace = parser.parse_args([''])  # empty jobStore attribute
    namespace.jobStore = None  # jobStore attribute will be updated per-batch
    return namespace
