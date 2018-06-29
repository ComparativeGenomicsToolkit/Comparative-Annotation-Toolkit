# Copyright 2006-2012 Mark Diekhans
"""
Wrapper for pipeline.py. Provides a simple interface for compli'cat'ed unix-style process pipes.
"""
import sys
import pipeline
import subprocess


def call_proc(cmd, keepLastNewLine=False):
    """call a process and return stdout, exception with stderr in message.
    The  cmd is either a list of command and arguments, or pipeline, specified by
    a list of lists of commands and arguments."""
    stdout = pipeline.DataReader()
    cmd = getDockerCommand('cat',cmd)
    pl = pipeline.Procline(cmd, stdin="/dev/null", stdout=stdout, stderr='/dev/stderr')
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
    cmd = getDockerCommand('cat',cmd)
    pl = pipeline.Procline(cmd, stdin=stdin, stdout=stdout, stderr=stderr)
    pl.wait()


def run_proc_code(cmd, stdin="/dev/null", stdout=None, stderr=None):
    """run a process, with I/O redirection to specified file paths or open
    file objects. None specifies inheriting open file.  Return exit code rather
    than raising exception"""
    try:
        cmd = getDockerCommand('cat',cmd)
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
    command = getDockerCommand('cat',cmd)
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

import os



def mrca_path(path1, path2):

    """

    Gives the Most Recent Common Ancestor directory that contains both paths.

    

    >>> mrca_path('/usr/lib/python2.7', '/usr/bin/python')

    '/usr/'

    >> mrca_path('/usr/', '/usr/')

    '/usr/'

    """

    while True:

        if path1 == path2:

            return path1 + os.path.sep

        elif len(path1) > len(path2):

            path1 = os.path.dirname(path1)

        else:

            path2 = os.path.dirname(path2)

    raise RuntimeError('something weird happened: the two paths %s and %s had no common prefix.' % (path1, path2))



def getDockerCommand(image, cmd):

    """

    Takes in a command (as a list of arguments like ['halStats',

    'file']) and outputs another list of arguments that will run it in

    the given Docker container, relativizing paths and binding

    directories when necessary.



    image: the Docker image to use, e.g. 'quay.io/comparative-genomics-toolkit/cactus:latest'

    cmd: list of arguments

    """

    dockerPreamble = ['docker', 'run', '-i', '--rm']

    # Find work_dir (MRCA of all provided files)

    work_dir = None

    for i, arg in enumerate(cmd):

        if os.path.isfile(arg):
            arg = os.path.abspath(arg)
            cmd[i] = arg
            if work_dir is None:

                work_dir = os.path.dirname(arg)

            else:

                work_dir = mrca_path(os.path.dirname(arg), work_dir)

    # Relativize all paths.

    if work_dir is not None:
        work_dir = os.path.abspath(work_dir)

        cmd = [arg.replace(work_dir, '/data/') for arg in cmd]

        dockerPreamble += ['-v', work_dir + ':/data']

    return dockerPreamble + [image] + cmd



