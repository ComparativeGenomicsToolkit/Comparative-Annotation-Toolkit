"""
Provides a simple interface between Toil and Luigi.
"""
import os
import luigi
import shutil
import bio
import fileOps
from toil.job import Job


class ToilOptionsMixin(object):
    """
    Add to a luigi WrapperTask that will be a entry point so that toil options can propagate downwards.
    """
    workDir = luigi.Parameter(default=None, significant=False)  # set up before use
    batchSystem = luigi.Parameter(default='singleMachine', significant=False)
    maxCores = luigi.IntParameter(default=16, significant=False)
    logLevel = luigi.Parameter(default='WARNING', significant=False)  # this is passed to toil
    cleanWorkDir = luigi.Parameter(default='onSuccess', significant=False)  # debugging option
    parasolCommand = luigi.Parameter(default=None, significant=False)


class ToilTask(luigi.Task, ToilOptionsMixin):
    """
    Task for launching toil pipelines from within luigi.
    """
    def prepare_toil_options(self, work_dir):
        """
        Prepares a Namespace object for Toil which has all defaults, overridden as specified
        Will see if the jobStore path exists, and if it does, assume that we need to add --restart
        :param work_dir: Parent directory where toil work will be done. jobStore will be placed inside. Will be used
        to fill in the workDir class variable.
        :return: Namespace
        """
        job_store = os.path.join(work_dir, 'jobStore')
        fileOps.ensure_file_dir(job_store)
        toil_args = get_toil_defaults()
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


def get_toil_defaults():
    """
    Extracts the default toil options as a dictionary, setting jobStore to None
    :return: dict
    """
    parser = Job.Runner.getDefaultArgumentParser()
    namespace = parser.parse_args([''])  # empty jobStore attribute
    namespace.jobStore = None  # jobStore attribute will be updated per-batch
    return namespace


###
# Helper functions for luigi-toil pipelines
###


def load_fasta_from_filestore(job, fasta_file_id, fasta_gdx_file_id, fasta_flat_file_id, prefix='genome', upper=False):
    """
    Convenience function that will load a fasta from the fileStore and return the local path to it. This works with
    the pyfasta module to load all of the required files.
    :param job: current job.
    :param fasta_file_id: fileStore file ID for the fasta
    :param fasta_gdx_file_id: fileStore file ID for the index (gdx)
    :param fasta_flat_file_id: fileStore file ID for the flat file (sentinel marker)
    :param prefix: local file path prefix
    :param upper: force all entries to upper case
    :return: open pyfasta Fasta record pointing to the file.
    """
    fasta_local_path = '{}.fasta'.format(prefix)
    job.fileStore.readGlobalFile(fasta_file_id, fasta_local_path)
    job.fileStore.readGlobalFile(fasta_gdx_file_id, '{}.fasta.gdx'.format(prefix))
    job.fileStore.readGlobalFile(fasta_flat_file_id, '{}.fasta.flat'.format(prefix))
    return bio.get_sequence_dict(fasta_local_path, upper=upper)


def write_fasta_to_filestore(toil, fasta_local_path):
    """
    Convenience function that loads a fasta and its associated gdx/flat file into the fileStore.
    Assumes that the paths are consistent with the requirements (i.e. $path.gdx and $path.flat)
    :param toil: Toil context manager
    :param fasta_local_path: Path to local fasta to load.
    :return: List of fileStore IDs for fasta, fasta_gdx, fasta_flat
    """
    fasta_file_id = toil.importFile('file:///' + fasta_local_path)
    gdx_file_id = toil.importFile('file:///' + fasta_local_path + '.gdx')
    flat_file_id = toil.importFile('file:///' + fasta_local_path + '.flat')
    return fasta_file_id, gdx_file_id, flat_file_id
