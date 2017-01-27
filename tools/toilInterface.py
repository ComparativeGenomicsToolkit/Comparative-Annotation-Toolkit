"""
Helper functions for toil-luigi interfacing
"""
import bio

###
# Helper functions for luigi-toil pipelines
###


def load_fasta_from_filestore(job, fasta_file_ids, prefix='genome', upper=False):
    """
    Convenience function that will load a fasta from the fileStore and return the local path to it. This works with
    the pyfasta module to load all of the required files.
    :param job: current job.
    :param fasta_file_ids: list of fileStore file ID for the fasta, gdx, and flat file.
    :param prefix: local file path prefix
    :param upper: force all entries to upper case
    :return: open pyfasta Fasta record pointing to the file.
    """
    fasta_local_path = '{}.fasta'.format(prefix)
    fasta_file_id, gdx_file_id, flat_file_id = fasta_file_ids
    job.fileStore.readGlobalFile(fasta_file_id, fasta_local_path)
    job.fileStore.readGlobalFile(gdx_file_id, '{}.fasta.gdx'.format(prefix))
    job.fileStore.readGlobalFile(flat_file_id, '{}.fasta.flat'.format(prefix))
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
