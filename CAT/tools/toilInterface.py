"""
Provides a simple interface between Toil and Luigi.
"""
import os
import tempfile
import luigi
import fileOps
from toil.job import Job


class ToilTask(luigi.Task):
    """
    Task for launching toil pipelines from within luigi.
    """
    workDir = luigi.Parameter(default=tempfile.gettempdir(), significant=False)
    batchSystem = luigi.Parameter(default='singleMachine', significant=False)
    maxCores = luigi.IntParameter(default=16, significant=False)
    logLevel = luigi.Parameter(default='WARNING', significant=False)

    def prepare_toil_options(self, job_store=None):
        """
        Prepares a Namespace object for Toil which has all defaults, overridden as specified
        Will see if the jobStore path exists, and if it does, assume that we need to add --restart
        :param job_store: path to jobStore. Will use default if not set.
        :return: Namespace
        """
        if job_store is None:
            job_store = os.path.abspath(os.path.join(self.workDir, tempfile.gettempdir()))
        else:
            job_store = os.path.abspath(job_store)
        toil_args = get_toil_defaults()
        toil_args.__dict__.update(vars(self))
        toil_args.jobStore = job_store
        if os.path.exists(job_store):
            toil_args.restart = True
        else:
            fileOps.ensure_file_dir(job_store)
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
