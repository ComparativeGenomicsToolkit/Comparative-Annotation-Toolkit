"""
File operations.

Original Author: Mark Diekhans
Modified: Ian Fiddes
"""

import os
import errno
import socket
import shutil
import gzip
import string
import random
import tempfile
import hashlib
import six


class TemporaryFilePath(object):
    """
    Generates a path pointing to a temporary file. Context manager wrapper for get_tmp_file. Deletes the file on exit.
    """
    def __init__(self, prefix=None, suffix="tmp", tmp_dir=None):
        self.path = get_tmp_file(prefix=prefix, suffix=suffix, tmp_dir=tmp_dir)

    def __enter__(self):
        return self.path

    def __exit__(self, type, value, traceback):
        try:
            os.remove(self.path)
        except OSError:
            pass


class TemporaryDirectoryPath(object):
    """
    Generates a path pointing to a temporary directory. Context manager wrapper for get_tmp_file,
    except creates a directory out of the path. Deletes the directory and all of its contents on exit.
    """
    def __init__(self, prefix=None, suffix="tmp", tmp_dir=None):
        self.path = get_tmp_file(prefix=prefix, suffix=suffix, tmp_dir=tmp_dir)
        ensure_dir(self.path)

    def __enter__(self):
        return self.path

    def __exit__(self, type, value, traceback):
        try:
            shutil.rmtree(self.path)
        except OSError:
            pass


def dir_is_writeable(d):
    """
    Is directory d writeable?
    :param d: directory
    :return: boolean
    """
    return os.access(d, os.W_OK | os.X_OK)


def ensure_dir(d):
    """
    Ensure that a directory exists, creating it (and parents) if needed.
    :param d: directory path to create.
    """
    try:
        os.makedirs(d)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(d):
            pass
        elif len(d) == 0:
            pass
        else:
            raise RuntimeError('Unable to create directory {}'.format(d))
    if not dir_is_writeable(d):
        raise RuntimeError('{} is not writeable.'.format(d))


def ensure_file_dir(file_path):
    """
    Ensure that the parent directory of a file path exists, creating it as needed, and making sure it is writeable.
    :param file_path: Path of file to ensure a parent directory of.
    """
    d = os.path.dirname(file_path)
    if d != '':
        ensure_dir(d)


def opengz(file, mode="r"):
    """
    Transparently open a gzipped or non-gzipped file for reading or writing. 
    :param file: Path of file to open for writing.
    :param mode: Same mode options as python's default open.
    :return: A open file handle.
    """
    assert mode in ['r', 'rb', 'a', 'ab', 'w', 'wb']
    if mode == 'wb' or (mode == 'w' and file.endswith('.gz')):
        return gzip.open(file, 'wb')
    elif mode == 'ab' or (mode == 'a' and file.endswith('.gz')):
        return gzip.open(file, 'ab')
    elif mode == 'w':
        return open(file, 'w')
    f = open(file, 'rb')
    if f.read(2) == '\x1f\x8b':
        f.seek(0)
        return gzip.GzipFile(fileobj=f, mode=mode)
    else:
        f.close()
        return open(file, mode)


def iter_lines(fspec, skip_lines=0, sep='\t'):
    """generator over lines in file, dropping newlines.  If fspec is a string,
    open the file and close at end. Otherwise it is file-like object and will
    not be closed.
    :param fspec: A file path or file handle.
    :param skip_lines: A integer of the number of lines to skip from the start of the file
    :param sep: Character used to separate columns in the file. If set to None, will not split the line.
    :return: Iterator of lines"""
    fh = _resolve_fspec(fspec, 'r')
    try:
        _ = [next(fh) for _ in range(skip_lines)]
        for line in fh:
            if sep is not None:
                yield line.rstrip().split(sep)
            else:
                yield line.rstrip()
    finally:
        if isinstance(fspec, six.string_types):
            fh.close()


def get_tmp_file(prefix=None, suffix="tmp", tmp_dir=None):
    """
    Returns the path to a temporary file. This file is guaranteed to not exist.
    :param prefix: Prefix to add to file path.
    :param suffix: Suffix to add to file path.
    :param tmp_dir: Directory to use. If None, will attempt to make use of system variables to find a path.
    :return: A file path.
    """
    if tmp_dir is None:
        tmp_dir = tempfile.gettempdir()
    if prefix is None:
        base_path = os.path.join(tmp_dir, '.'.join([socket.gethostname(), str(os.getpid())]))
    else:
        base_path = os.path.join(tmp_dir, '.'.join([prefix, socket.gethostname(), str(os.getpid())]))
    while True:
        rand = ''.join([random.choice(string.digits) for _ in range(10)])
        path = '.'.join([base_path, rand, suffix])
        if not os.path.exists(path):
            return path


def get_tmp_toil_file(prefix=None, suffix="tmp"):
    """
    Returns the path to a temporary file. This is a convenience wrapper for get_tmp_file that sets tmp_dir to
    os.getcwd().
    This is useful because of how toil caching works.
    It also returns the absolute path.

    :param prefix: Prefix to add to file path.
    :param suffix: Suffix to add to file path.
    :return: A file path.
    """
    return os.path.abspath(get_tmp_file(prefix=prefix, suffix=suffix, tmp_dir=os.getcwd()))


def atomic_install(tmp_path, final_path):
    """
    Atomically install a file from tmp_path to final_path. Handles crossing file system boundaries.
    :param tmp_path: Path of parent (temporary) file.
    :param final_path: Destination path.
    """
    ensure_file_dir(final_path)
    try:
        os.rename(tmp_path, final_path)
    except OSError:
        tmp = get_tmp_file(tmp_dir=os.path.dirname(final_path))
        shutil.copy(tmp_path, tmp)
        os.rename(tmp, final_path)
        os.remove(tmp_path)


def touch(file_path):
    """
    Creates a blank file at file path, ensuring it exists
    :param file_path: string to file
    :return: None
    """
    ensure_file_dir(file_path)
    with open(file_path, 'a'):
        os.utime(file_path, None)


def print_row(fspec, line, sep='\t'):
    """
    Convenience function that writes a delimited line to fspec (file handle or file)
    :param fspec: A open file handle or file path
    :param line: One or more things to write. Must be convertible to strings.
    :param sep: separator to use
    """
    fh = _resolve_fspec(fspec, 'w')
    fh.write(sep.join(map(str, line)) + '\n')


def print_rows(fspec, item_iter, sep='\t'):
    """
    Convenience function that writes a iterable of lines to fspec (file handle or file)
    :param fspec: A open file handle or file path
    :param item_iter: One or more things to write. Must be convertible to strings.
    :param sep: separator to use
    """
    fh = _resolve_fspec(fspec, 'w')
    for line in item_iter:
        print_row(fh, line, sep)


def print_iterable(fspec, item_iter):
    """
    Convenience function that simply writes an iterable of lines to fspec (file handle or file)
    :param fspec: A open file handle or file path
    :param item_iter: One or more things to write. Assumed to be fully formatted strings with newlines
    """
    fh = _resolve_fspec(fspec, 'w')
    for line in item_iter:
        fh.write(line)


def _resolve_fspec(fspec, mode='r'):
    """
    Determine if this is a file or a handle, passing a file name to opengz()
    :param fspec: A open file handle or file path
    :return: a open file handle
    """
    if isinstance(fspec, six.string_types):
        return opengz(fspec, mode)
    else:
        return fspec


def hashfile(fspec, hasher=hashlib.sha256, blocksize=65536, num_characters=12):
    """
    Calculates a SHA256 hash of a file.
    :param fspec: path or handle
    :param hasher: hashing function to use
    :param blocksize: size of file blocks to work on
    :param num_characters: Number of characters to use of hex digest
    :return: integer
    """
    fh = _resolve_fspec(fspec)
    buf = fh.read(blocksize)
    hasher = hasher()  # instantiate this hashing instance
    while len(buf) > 0:
        hasher.update(buf.encode('utf-8'))
        buf = fh.read(blocksize)
    return int(hasher.hexdigest(), 16) % 10 ** num_characters
