"""
Helper functions for toil-luigi interfacing
"""
from . import bio
import math
import argparse
from toil.lib.humanize import human2bytes
try:
    from toil.fileStores import FileID
except ImportError:
    from toil.fileStore import FileID

###
# Helper functions for luigi-toil pipelines
###


def load_fasta_from_filestore(job, fasta_file_ids, prefix='genome', upper=False):
    """
    Convenience function that will load a fasta from the fileStore and return the local path to it. This works with
    the pyfaidx module to load all of the required files.
    :param job: current job.
    :param fasta_file_ids: list of fileStore file ID for the fasta and fasta index files.
    :param prefix: local file path prefix
    :param upper: force all entries to upper case
    :return: open pyfaidx Fasta record pointing to the file.
    """
    fasta_local_path = '{}.fasta'.format(prefix)
    fasta_file_id, fai_file_id = fasta_file_ids
    job.fileStore.readGlobalFile(fasta_file_id, fasta_local_path)
    job.fileStore.readGlobalFile(fai_file_id, '{}.fasta.fai'.format(prefix))
    return bio.get_sequence_dict(fasta_local_path, upper=upper)


def write_fasta_to_filestore(toil, fasta_local_path):
    """
    Convenience function that loads a fasta and its associated gdx/flat file into the fileStore.
    Assumes that the paths are consistent with the requirements (i.e. $path.fai)
    :param toil: Toil context manager
    :param fasta_local_path: Path to local fasta to load.
    :return: List of fileStore IDs for fasta, fasta_fai
    """
    fasta_file_id = FileID.forPath(toil.importFile('file:///' + fasta_local_path), fasta_local_path)
    fai_file_id = FileID.forPath(toil.importFile('file:///' + fasta_local_path + '.fai'), fasta_local_path + '.fai')
    return fasta_file_id, fai_file_id


def find_total_disk_usage(input_file_ids, buffer='2G', round='2G'):
    """
    Takes a input_file_id namespace or dict or list and finds all members that are FileID objects,
    and finds their sizes.
    Based on buffer and round, returns a integer value of disk usage in bytes to pass to a toil job.
    :param input_file_ids: A namespace object with an arbitrary nesting of possible file ID values
    :param buffer: Additional space buffer requested. Human readable parsed by human2bytes
    :param round: amount to round up. Human readable parsed by human2bytes
    :return: integer
    """
    def roundup(x, base):
        return int(math.ceil(x / float(base))) * base

    def descend_object(obj):
        if isinstance(obj, dict):
            for item in list(obj.values()):
                for v in descend_object(item):
                    yield v
        elif isinstance(obj, list):
            for item in obj:
                for v in descend_object(item):
                    yield v
        elif isinstance(obj, argparse.Namespace):
            for item in list(obj.__dict__.values()):
                for v in descend_object(item):
                    yield v
        elif isinstance(obj, FileID):
            yield obj

    tot = sum([x.size for x in descend_object(input_file_ids)])
    return roundup(tot, human2bytes(round)) + human2bytes(buffer)
