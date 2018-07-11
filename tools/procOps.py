# Copyright 2006-2012 Mark Diekhans
"""
Wrapper for pipeline.py. Provides a simple interface for compli'cat'ed unix-style process pipes.
"""
import os
import pipeline
import sys
import subprocess

def cmdLists(cmd):
    """creates dockers commands from lists or a list of lists
    """
    if isinstance(cmd[0],list):
        docList = []
        for e in cmd:
            docList.append(getDockerCommand('cat',e))
        return docList
    else:
        return getDockerCommand('cat',cmd)


def call_proc(cmd, keepLastNewLine=False):
    """call a process and return stdout, exception with stderr in message.
    The  cmd is either a list of command and arguments, or pipeline, specified by
    a list of lists of commands and arguments."""
    stdout = pipeline.DataReader()
    cmd = cmdLists(cmd)
    pl = pipeline.Procline(cmd, stdin="/dev/null", stdout=stdout)
    pl.wait()
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
    dockerPreamble = ['docker', 'run', '-i', '--rm']
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
