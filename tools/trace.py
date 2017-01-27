# Copyright 2006-2012 Mark Diekhans
"""Debug tracing to a file"""
# ideas from:
#  http://www.dalkescientific.com/writings/diary/archive/2005/04/20/tracing_python_code.html

import sys, os, types, linecache, threading

# used to detect traces that are currently open to prevent closing
_activeTraceFds = set()

def getActiveTraceFds():
    "return snapshot of currently active traces"
    return frozenset(_activeTraceFds)

class Trace(object):
    """Trace object, associate with an open trace file.  File is flushed after
    each write to debug blocking"""

    def __init__(self, traceFile, ignoreMods=None, inclThread=False, inclPid=False, callIndent=True):
        "open log file. ignoreMods is a sequence of module objects or names to skip"
        self.fh = open(traceFile, "w")
        _activeTraceFds.add(self.fh.fileno())
        self.inclThread = inclThread
        self.inclPid = inclPid
        self.callIndent = callIndent
        self.depth = 0
        self.ignoreMods = set()
        if ignoreMods is not None:
            for m in ignoreMods:
                if type(m) == types.ModuleType:
                    m = m.__name__
                self.ignoreMods.add(m)

    def enable(self):
        """enable logging on all threads."""
        assert(self.fh is not None)
        sys.settrace(self.__callback)
        threading.settrace(self.__callback)

    def disable(self):
        """disable logging on all threads."""
        sys.settrace(None)
        threading.settrace(None)

    def close(self):
        "disable and close log file"
        if self.fh is not None:
            self.disable()
            _activeTraceFds.remove(self.fh.fileno())
            self.fh.close()

    def log(self, *args):
        """log arguments as a message followed by a newline"""
        # build string and output, as this minimize interleaved messages
        # discard output on I/O errors
        try:
            msg = []
            if self.inclPid:
                msg.append(str(os.getpid()))
                msg.append(": ")
            if self.inclThread:
                # can only include id, getting name will cause deadlock
                msg.append(str(threading._get_ident()))
                msg.append(": ")
            for a in args:
                msg.append(str(a))
            msg.append("\n")
            self.fh.write("".join(msg))
            self.fh.flush()
        except IOError as ex:
            pass

    __indentStrs = {} # cache of spaces for identation, indexed by depth
    def __getIndent(self):
        "get indentation string"
        if not self.callIndent:
            return ""
        i = Trace.__indentStrs.get(self.depth)
        if i is None:
            i = Trace.__indentStrs[self.depth] = "".ljust(4*self.depth)
        return i

    def __logLine(self, frame, event):
        "log a code line"
        lineno = frame.f_lineno
        fname = frame.f_globals["__file__"]
        if (fname.endswith(".pyc") or fname.endswith(".pyo")):
            fname = fname[:-1]
        name = frame.f_globals["__name__"]
        line = linecache.getline(fname, lineno)
        self.log(name, ":", lineno, self.__getIndent(), line.rstrip())

    __logEvents = frozenset(["call", "line"])
    def __callback(self, frame, event, arg):
        "trace event callback"
        if frame.f_globals["__name__"] not in self.ignoreMods:
            if event == "call":
                self.depth += 1
            elif event == "return":
                self.depth -= 1
                if self.depth < 0:
                    self.depth = 0
            
            if event in Trace.__logEvents:
                self.__logLine(frame, event)
        return self.__callback

__all__ = (getActiveTraceFds.__name__, Trace.__name__)
