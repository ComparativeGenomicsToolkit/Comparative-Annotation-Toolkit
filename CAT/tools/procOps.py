# Copyright 2006-2012 Mark Diekhans
"""
Wrapper for pipeline.py. Provides a simple interface for complicated unix-style process pipes.
"""
import sys
import pipeline
import subprocess


def call_proc(cmd, keepLastNewLine=False):
    """call a process and return stdout, exception with stderr in message.
    The  cmd is either a list of command and arguments, or pipeline, specified by
    a list of lists of commands and arguments."""
    stdout = pipeline.DataReader()
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
    pl = pipeline.Procline(cmd, stdin=stdin, stdout=stdout, stderr=stderr)
    pl.wait()


def run_proc_code(cmd, stdin="/dev/null", stdout=None, stderr=None):
    """run a process, with I/O redirection to specified file paths or open
    file objects. None specifies inheriting open file.  Return exit code rather
    than raising exception"""
    try:
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
    if stdin != None:
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
