# Copyright 2006-2012 Mark Diekhans
"""
Wrapper for pipeline.py. Provides a simple interface for compli'cat'ed unix-style process pipes.
"""
import os
import pipeline
import pipes
import sys
import subprocess
import logging
import time

logger = logging.getLogger(__name__)

def cmdLists(cmd):
    """
    creates docker or singularity command(s) from either a
    single command or a list of commands.
    """
    if os.environ.get('CAT_BINARY_MODE') == 'docker':
        if isinstance(cmd[0],list):
            docList = []
            for e in cmd:
                docList.append(getDockerCommand('quay.io/ucsc_cgl/cat',e))
            return docList
        else:
            return getDockerCommand('quay.io/ucsc_cgl/cat',cmd)
    elif os.environ.get('CAT_BINARY_MODE') == 'singularity':
        img = os.path.join(os.environ.get('SINGULARITY_PULLFOLDER'), 'cat.img')
        if isinstance(cmd[0], list):
            return list(map(lambda c: getSingularityCommand(img, c), cmd))
        else:
            return getSingularityCommand(img, cmd)
    else:
        return cmd

def call_proc(cmd, keepLastNewLine=False):
    """call a process and return stdout, exception with stderr in message.
    The  cmd is either a list of command and arguments, or pipeline, specified by
    a list of lists of commands and arguments."""
    stdout = pipeline.DataReader()
    cmd = cmdLists(cmd)
    logger.debug('About to run command: %s' % cmd)
    now = time.time()
    pl = pipeline.Procline(cmd, stdin="/dev/null", stdout=stdout)
    pl.wait()
    logger.debug('Command %s took %s seconds.' % (cmd, time.time() - now))
    out = stdout.get()
    if (not keepLastNewLine) and (len(out) > 0) and (out[-1] == "\n"):
        out = out[0:-1]
    return out


def call_proc_lines(cmd):
    """call a process and return stdout, split into a list of lines, exception
    with stderr in message."""
    out = call_proc(cmd)
    if len(out) == 0:
        return []  # split creates a list of one empty string from an empty string
    return out.split("\n")


def run_proc(cmd, stdin="/dev/null", stdout=None, stderr=None):
    """run a process, with I/O redirection to specified file paths or open
    file objects. None specifies inheriting open file."""
    cmd = cmdLists(cmd)
    pl = pipeline.Procline(cmd, stdin=stdin, stdout=stdout, stderr=stderr)
    pl.wait()


def run_proc_code(cmd, stdin="/dev/null", stdout=None, stderr=None):
    """run a process, with I/O redirection to specified file paths or open
    file objects. None specifies inheriting open file.  Return exit code rather
    than raising exception"""
    try:
        cmd = cmdLists(cmd)
        pl = pipeline.Procline(cmd, stdin=stdin, stdout=stdout, stderr=stderr)
        pl.wait()
    except pipeline.ProcException as ex:
        if ex.returncode is not None:
            return ex.returncode
        else:
            raise ex
    return 0


def popen_catch(command, stdin=None):
    """
    Runs a command and return standard out. TODO: use Mark's tools. I don't think he has this functionality.
    """
    command = cmdLists(command)
    if stdin is not None:
        process = subprocess.Popen(command,
                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=sys.stderr, bufsize=-1)
        output, nothing = process.communicate(stdin)
    else:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=sys.stderr, bufsize=-1)
        output, nothing = process.communicate()
    sts = process.wait()
    if sts != 0:
        raise RuntimeError("Command: %s with stdin string '%s' exited with non-zero status %i" % (command, stdin, sts))
    return output


def mrca_path(path1, path2):
    """
    Gives the Most Recent Common Ancestor directory that contains both paths.
    
    >>> mrca_path('/usr/lib/python2.7', '/usr/bin/python')
    '/usr'
    >>> mrca_path('/usr/', '/usr/')
    '/usr'
    """
    while True:
        if path1 == path2:
            return os.path.normpath(path1)
        elif len(path1) > len(path2):
            path1 = os.path.dirname(path1)
        else:
            path2 = os.path.dirname(path2)
    raise RuntimeError('something weird happened: the two paths %s and %s had no common prefix.' % (path1, path2))

def add_to_work_dirs(dirname, work_dirs):
    """
    >>> work_dirs = []
    >>> add_to_work_dirs('/tmp', work_dirs)
    >>> work_dirs
    ['/tmp']
    >>> add_to_work_dirs('/foo/bar/baz', work_dirs)
    >>> work_dirs
    ['/tmp', '/foo/bar/baz']
    >>> add_to_work_dirs('/foo/baz', work_dirs)
    >>> work_dirs
    ['/tmp', '/foo']
    """
    if not work_dirs:
        # Empty list
        work_dirs.append(dirname)
    else:
        # Go through all existing work_dirs and see if there's one
        # that can be added without traversing back to the root. (We
        # don't want to bind-mount the root dir, because that will
        # override the existing filesystem within the container.)
        for i, work_dir in enumerate(work_dirs):
            mrca = mrca_path(dirname, work_dir)
            if mrca == '/':
                # Avoid bind-mounting the root dir.
                if i == len(work_dirs) - 1:
                    # No mergeable directories.
                    work_dirs.append(dirname)
            else:
                # Mergable. Replace this entry with the MRCA.
                work_dirs[i] = mrca

def getDockerCommand(image, cmd):
    """
    Takes in a command (as a list of arguments like ['halStats',
    'file']) and outputs another list of arguments that will run it in
    the given Docker container, binding directories when necessary.

    image: the Docker image to use, e.g. 'quay.io/comparative-genomics-toolkit/cactus:latest'
    cmd: list of arguments
    """
    dockerPreamble = ['docker', 'run', '-i', '--rm', '-u', "%s:%s" % (os.getuid(), os.getgid())]
    work_dirs = []
    for i, arg in enumerate(cmd):
        if arg.startswith('-') and '=' in arg:
            # We assume this is -option=value syntax. Special-case
            # this to check if the value is a path.
            arg = arg.split('=')[1]
        dirname = os.path.dirname(arg)
        if os.path.exists(dirname):
            # The dirname exists, so we will try to mount it.
            arg = os.path.abspath(arg)
            if arg.startswith('/dev'):
                continue
            add_to_work_dirs(dirname, work_dirs)
    for work_dir in work_dirs:
        work_dir = os.path.abspath(work_dir)
        dockerPreamble += ['-v', work_dir + ':' + work_dir]
    return dockerPreamble + [image] + cmd

def getSingularityCommand(image, cmd):
    """
    Takes a command and turns it into a Singularity command
    that will run the original command inside a container.

    image: the Singularity image to use. For some reason,
        it is much faster to give Singularity a pre-built
        '.img' file than to give it a URL, even if it has
        already cached that URL.
    cmd: command to turn into a Singularity command,
        represented as a list of arguments
    """
    # singularity only allows mounting at existing mount
    # points by default, so rather than mounting the
    # directory for each file argument as in
    # getDockerCommand, we mount the entire root of the
    # outside file system in '/mnt' of the container, and
    # then prepend '/mnt' to all file paths in the command.
    singularity_cmd = ['singularity', 'exec', '-B', '/:/mnt', image]
    for arg in cmd:
        if arg.startswith('-') and len(arg.split('=')) == 2:
            # We assume this is -option=value syntax. Special-case
            # this to check if the value is a path.
            option, value = arg.split('=')
            singularified_arg = '='.join([option, singularify_arg(value)])
        else:
            singularified_arg = singularify_arg(arg)

        singularity_cmd.append(singularified_arg)

    return singularity_cmd

def singularify_arg(arg, singularity_mount_point='/mnt'):
    """
    Check to see if 'arg' is a path; if it is, modify it to
    be accessible from inside the singularity container.

    arg: String containing an argument to singularify
    singularity_mount_point: the place where the outside
        root is mounted in the singularity container

    Returns: arg if arg is not a path, or a container-
        accessible version of arg if it is a path
    """
    # first, check to see if this argument is a file,
    # or could potentially be a file in the future, by
    # seeing if its dirname exists
    if os.path.exists(os.path.dirname(arg)):
        # prepend /mnt to the path because we mounted
        # '/' of the outside to '/mnt' of the container
        # (os.path.join() cannot be used to prepend
        # like this)
        arg = '/mnt/' + os.path.abspath(arg)

    return arg

