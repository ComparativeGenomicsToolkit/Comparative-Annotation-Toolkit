"""
Provides a simple interface between Toil and Luigi.
"""
import os
import tempfile
import luigi
import bio
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
    cleanWorkDir = luigi.Parameter(default=None, significant=False)

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


def load_fasta_from_filestore(job, fasta_file_id, fasta_gdx_file_id, fasta_flat_file_id, prefix='genome',
                             upper=False):
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
