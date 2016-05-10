# Copyright 2006-2012 Mark Diekhans
"""
Wrapper for pipeline.py. Provides a simple interface for complicated unix-style process pipes.
"""

import tools.pipeline


def call_proc(cmd, keepLastNewLine=False):
    """call a process and return stdout, exception with stderr in message.
    The  cmd is either a list of command and arguments, or pipeline, specified by
    a list of lists of commands and arguments."""
    stdout = tools.pipeline.DataReader()
    pl = tools.pipeline.Procline(cmd, stdin="/dev/null", stdout=stdout)
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
    pl = tools.pipeline.Procline(cmd, stdin=stdin, stdout=stdout, stderr=stderr)
    pl.wait()


def run_proc_code(cmd, stdin="/dev/null", stdout=None, stderr=None):
    """run a process, with I/O redirection to specified file paths or open
    file objects. None specifies inheriting open file.  Return exit code rather
    than raising exception"""
    try:
        pl = tools.pipeline.Procline(cmd, stdin=stdin, stdout=stdout, stderr=stderr)
        pl.wait()
    except tools.pipeline.ProcException as ex:
        if ex.returncode is not None:
            return ex.returncode
        else:
            raise ex
    return 0
