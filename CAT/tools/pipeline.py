# Copyright 2006-2012 Mark Diekhans
"""
Process pipelines constructed as a DAG.
"""
import os, sys, fcntl, signal, errno, threading, traceback, pickle, time
from tools import fifo, trace, strOps, PycbioException

try:
    MAXFD = os.sysconf("SC_OPEN_MAX")
except OSError:
    MAXFD = 256


def _getSigName(num):
    "get name for a signal number"
    # find name in signal namespace
    for key in vars(signal):
        if (getattr(signal, key) == num) and key.startswith("SIG") and (key.find("_") < 0):
            return key
    return "signal"+str(num)

def _setPgid(pid, pgid):
    """set pgid of a process, ignored exception caused by race condition
    that occurs if already set by parent or child has already existed"""
    # Should just ignore on EACCES, as to handle race condition with parent
    # and child.  However some Linux kernels (seen in 2.6.18-53) report ESRCH
    # or EPERM.  To handle this is a straight-forward way, just check that the
    # change has been made.  However, in some cases the change didn't take,
    # retrying seems to make the problem go away.
    for i in xrange(0,5):
        try:
            os.setpgid(pid, pgid)
            return
        except OSError:
            if os.getpgid(pid) == pgid:
                return
            time.sleep(0.25) # sleep for retry
    # last try, let it return an error
    os.setpgid(pid, pgid)

# FIXME: why not use pipes.quote?
def _quoteStr(a):
    "return string  with quotes if it contains white space"
    a = str(a).replace('"', '\\"')
    if strOps.hasSpaces(a):
        a = '"' + a + '"'
    return a

class ProcException(PycbioException):
    "Process error exception.  A None returncode indicates a exec failure."
    def __init__(self, procDesc, returncode=None, stderr=None, cause=None):
        self.returncode = returncode
        self.stderr = stderr
        if returncode is None:
            msg = "exec failed"
        elif (returncode < 0):
            msg = "process signaled: " + _getSigName(-returncode)
        else:
            msg = "process exited " + str(returncode)
        if procDesc is not None:
            msg += ": " + procDesc
        if (stderr is not None) and (len(stderr) != 0):
            msg += ":\n" + stderr
        PycbioException.__init__(self, msg, cause=cause)

class ProcDagException(PycbioException):
    "Exception not associate with process execution"
    def __init__(self, msg, cause=None):
        PycbioException.__init__(self, msg, cause)

def nonBlockClear(fd):
    "clear the non-blocking flag on a fd"
    flags = fcntl.fcntl(fd, fcntl.F_GETFL)
    fcntl.fcntl(fd, fcntl.F_SETFL, flags&~os.O_NONBLOCK)

class _StatusPipe(object):
    """Used to communicate from child to parent.  Child is close-on-exec,
    so this can be used to get the status of an exec."""
    __slots__ = ("rfd", "wfd")

    def __init__(self):
        self.rfd, self.wfd = os.pipe()
        flags = fcntl.fcntl(self.wfd, fcntl.F_GETFD)
        fcntl.fcntl(self.wfd, fcntl.F_SETFD, flags|fcntl.FD_CLOEXEC)

    def postForkParent(self):
        "post fork handling in parent"
        os.close(self.wfd)
        self.wfd = None

    def postForkChild(self):
        "post fork handling in child"
        os.close(self.rfd)
        self.rfd = None

    def sendExcept(self, ex):
        "send an exception to parent and exit child"
        os.write(self.wfd, pickle.dumps(ex))
        os._exit(255)

    def sendOk(self, ex):
        "send True to parent and continue"
        os.write(self.wfd, pickle.dumps(True))

    def recvStatus(self):
        """read status from child, return exception if received, otherwise
        None or True"""
        # FIXME add read loop, or read through pickle??
        data = os.read(self.rfd, 1024*1024)
        os.close(self.rfd)
        if len(data) > 0:
            return pickle.loads(data)
        else:
            return None 

class PInOut(object):
    """base class for PIn and POut"""
    def __init__(self, dev, argPrefix=None):
        self.dev = dev
        self.argPrefix = argPrefix
        self.proc = None
        self.named = None
        dev._addPio(self)

    def __radd__(self, argPrefix):
        "string concatiation operator that sets the argPrefix"
        assert(self.argPrefix is None)
        self.argPrefix = argPrefix
        return self

    def assocByPath(self, proc):
        "associate Dev with a Proc that will access it by path"
        assert(self.proc is None)
        self.proc = proc
        self.named = True

    def assocByFd(self, proc):
        "associate Dev with a Proc that will access it by file descriptor"
        assert(self.proc is None)
        self.proc = proc
        self.named = False

    def isPipe(self):
        "is this associate with a Pipe"
        return isinstance(self.dev, Pipe)

    def getConnectedProc(self):
        "If this is a Pipe to a proc, return proc, otherwise None"
        if isinstance(self.dev, Pipe):
            return self.dev.pout.proc if isinstance(self, PIn) else self.dev.pin.proc
        else:
            return None

    def getFd(self):
        "get file descriptor for this object"
        assert(not self.named)
        return self.dev.getFd(self)
        
    def getFh(self):
        "get file object for this object, or error if not supported by Dev"
        assert(not self.named)
        return self.dev.getFh(self)
        
    def getPath(self):
        "get path for this object"
        return self.dev.getPath(self)

    def getArg(self):
        "get path with optional argument prefix"
        if self.argPrefix is None:
            return self.getPath()
        else:
            return self.argPrefix + self.getPath()

    def close(self):
        "terminate association with device"
        self.dev.close(self)
        
    def __str__(self):
        """return input file argument"""
        if not self.named:
            return str(self.dev.getFd(self))
        elif self.argPrefix is None:
            return self.dev.getPath(self)
        else:
            return self.argPrefix + self.dev.getPath(self)

    @staticmethod
    def pIsPipe(obj):
        "test if obj is a PInOut Pipe object"
        return isinstance(obj, PInOut) and obj.isPipe()

    @staticmethod
    def pHasProc(obj):
        "test if obj is a PInOut object with a proc"
        return isinstance(obj, PInOut) and (obj.proc is not None)

    @staticmethod
    def pHasOtherProc(obj):
        "test if obj is a PInOut object with a proc on the other end"
        if isinstance(obj, PIn):
            return PInOut.pHasProc(obj.dev.pout)
        elif isinstance(obj, POut):
            return PInOut.pHasProc(obj.dev.pin)
        else:
            return False

class PIn(PInOut):
    """Process input object that links Dev object as input to a process,
    either as stdin or as a command line argument.  That is, it's output
    from the Dev and input to the Proc.

    The argPrefix attribute is used in construction arguments when an
    option-equals is prepended to the file (--in=fname).  The argPrefix can be
    specified as an option to the constructor, or in a string concatination
    ("--in="+PIn(d))."""
    def __init__(self, dev, argPrefix=None):
        PInOut.__init__(self, dev, argPrefix)

class POut(PInOut):
    """Process output object that links Dev object as output from a process,
    either as stdout/stderr or as a command line argument.  That is, it's input
    to the Dev and output from the Proc.

    The argPrefix attribute is used in construction arguments when an
    option-equals is prepended to the file (--out=fname).  The argPrefix can be
    specified as an option to the constructor, or in a string concatenation
    ("--out="+POut(d)).

    If append is True, file is opened with append access, if approriate.
    """
    def __init__(self, dev, argPrefix=None, append=False):
        PInOut.__init__(self, dev, argPrefix)
        self.append = append  # FIXME: not implemented

class Dev(object):
    """Base class for objects specifiying process input or output.  Usually
    implemented as pipes or named pipes, they provide a way of hide details of
    setting up interprocess communication.  Objects of this type are linked to
    processes by PIn or POut objects. """

    def __init__(self):
        self.pin = None   # dev output/process input
        self.pout = None  # dev input/process output

    def _addPio(self, pio):
        "add a PInOut object"
        if isinstance(pio, PIn):
            if (self.pin is not None) and (self.pin != pio):
                raise Exception("PIn object already associated with Dev")
            self.pin = pio
        elif isinstance(pio, POut):
            if (self.pout is not None) and (self.pout != pio):
                raise Exception("POut object already associated with Dev")
            self.pout = pio
        else:
            raise Exception("attempt to add object of unknown class to a Dev: " + str(pio.__class__))

    def needNamed(self):
        """does this device need a named pipe?"""
        return (((self.pin is not None) and (self.pin.named))
                or ((self.pout is not None) and (self.pout.named)))

    def preFork(self):
        """pre-fork setup."""
        pass

    def postExecParent(self):
        "called do any post-exec handling in the parent"
        pass
        
    def close(self, pio):
        """remove association of process with device; PInOut association remains
        for debugging purposes"""
        pass

    def finish(self):
        "called in parent when processing is complete"
        pass

    def getFd(self, pio):
        "get file descriptor for given PInOut object"
        raise AttributeError("getFd not implemented")

    def getFh(self, pio):
        "get file object for given PInOut object, or error if not supported"
        raise AttributeError("getFh not supported for this Dev: " + str(self.__class__))
        
    def getPath(self, pio):
        "get path for given PInOut object"
        raise AttributeError("getPath not implemented")
        
class DataReader(Dev):
    """Object to read data from process into memory via a pipe."""

    # FIXME make sure it can handled binary data
    def __init__(self):
        Dev.__init__(self)
        self.fifo = None
        self.data = []
        self.thread = None

    def __del__(self):
        "finalizer"
        if self.thread is not None:
            try:
                self.thread.join()
            except: pass
        if self.fifo is not None:
            try:
                self.fifo.close()
            except: pass

    def __str__(self):
        return "[DataWriter]"

    def preFork(self):
        "pre-fork setup"
        if self.pin is not None:
            raise Exception(self.__class__.__name__ + " can't be used for process input")
        self.fifo = fifo.factory()

    def postExecParent(self):
        "called to do any post-exec handling in the parent"
        if not self.pout.named:
            self.fifo.wclose()
        self.thread = threading.Thread(target=self.__reader)
        self.thread.start()
        
    def finish(self):
        "called in parent when processing is complete"
        if self.fifo is not None:
            self.fifo.wclose()  # call before join thread so it sees EOF
            if self.thread is not None:
                self.thread.join()
                self.thread = None
            self.fifo.close()

    def __reader(self):
        "child read thread function"
        self.data.append(self.fifo.getRfh().read())

    def get(self):
        "return buffered data as a string"
        # FIXME: does this work for binary?
        return "".join(self.data)

    def getFd(self, pio):
        "get file descriptor for given PInOut object"
        assert(pio == self.pout)
        return self.fifo.wfd
        
    def getPath(self, pio):
        "get path for given PInOut object"
        assert(pio == self.pout)
        return self.fifo.wpath
        
class DataWriter(Dev):
    """Object to write data from memory to process via a pipe."""

    def __init__(self, data):
        Dev.__init__(self)
        self.fifo = None
        self.data = data
        self.thread = None

    def __del__(self):
        "finalizer"
        if self.thread is not None:
            try:
                self.thread.join()
            except: pass
        if self.fifo is not None:
            try:
                self.fifo.close()
            except: pass

    def __str__(self):
        return "[DataWriter]"

    def preFork(self):
        "pre-fork setup"
        if self.pout is not None:
            raise Exception(self.__class__.__name__ + " can't be used for process output")
        self.fifo = fifo.factory()

    def postExecParent(self):
        "called to do any post-exec handling in the parent"
        if not self.pin.named:
            self.fifo.rclose()
        self.thread = threading.Thread(target=self.__writer)
        self.thread.start()
        
    def finish(self):
        "called in parent when processing is complete"
        if self.thread is not None:
            self.thread.join()
            self.thread = None
        if self.fifo is not None:
            self.fifo.close()

    def __writer(self):
        "write thread function"
        try:
            self.fifo.getWfh().write(self.data)
            self.fifo.wclose()
        except IOError as ex:
            if ex.errno != errno.EPIPE:
                raise

    def getFd(self, pio):
        "get file descriptor for given PInOut object"
        assert(pio == self.pin)
        return self.fifo.rfd
        
    def getPath(self, pio):
        "get path for given PInOut object"
        assert(pio == self.pin)
        return self.fifo.rpath
        
class Pipe(Dev):
    """Interprocess communication between two Procs, either by named or
    anonymous pipes.  One end can also be attached to read/write
    to/from the process pipeline."""

    def __init__(self):
        Dev.__init__(self)
        self.fifo = None

    def __del__(self):
        "finalizer"
        if self.fifo is not None:
            try:
                self.fifo.close()
            except: pass

    def __str__(self):
        return "[Pipe]"

    def preFork(self):
        "pre-fork setup"
        self.fifo = fifo.factory()

    def postExecParent(self):
        "called to do any post-exec handling in the parent"
        # only close if side is associated with a pipe to a proc, keep
        # open if pipe is to use
        if (not self.pout.named) and (self.pout.proc is not None):
            self.fifo.wclose()
        if (not self.pin.named) and (self.pin.proc is not None):
            self.fifo.rclose()
        
    def close(self, pio):
        """remove association of process with device; PInOut association remains
        for debugging purposes"""
        if pio == self.pin:
            self.fifo.rclose()
        elif pio == self.pout:
            self.fifo.wclose()
        else:
            assert(False)

    def finish(self):
        "called in parent when processing is complete"
        if self.fifo is not None:
            self.fifo.close()

    def getFd(self, pio):
        "get file descriptor for given PInOut object"
        if pio == self.pin:
            return self.fifo.rfd
        else:
            return self.fifo.wfd
        
    def getFh(self, pio):
        "get file object for given PInOut object"
        if pio == self.pin:
            return self.fifo.getRfh()
        else:
            return self.fifo.getWfh()

    def getPath(self, pio):
        "get path for given PInOut object"
        if isinstance(pio, PIn):
            return self.fifo.rpath
        else:
            return self.fifo.wpath
        
class File(Dev):
    """A file path for input or output, used for specifying stdio associated
    with files. Proc wraps these around string arguments automatically"""

    def __init__(self, path, append=False):
        """constructor"""
        Dev.__init__(self)
        self.path = path
        self.append = append
        self.fd = None  # in child only

    def __str__(self):
        return self.path

    def getFd(self, pio):
        "get file descriptor for given PInOut object"
        if self.fd is None:
            if isinstance(pio, PIn):
                self.fd = os.open(self.path, os.O_RDONLY)
            elif self.append:
                self.fd = os.open(self.path, os.O_WRONLY|os.O_CREAT|os.O_APPEND, 0666)
            else:
                self.fd = os.open(self.path, os.O_WRONLY|os.O_CREAT|os.O_TRUNC, 0666)
        return self.fd
        
    def getPath(self, pio):
        "get path for given PInOut object"
        return self.path

class Proc(object):
    """A process, represented as a node in a DAG of Proc objects, connected by
    PInOut and Dev objects.  All processes in a ProcDag are part of the same
    process group.

    Process arguments can be PIn or POut object, or any other object, in which
    case they are converted by str() before exec

    If the stdin/out/err arguments can have the following values:
       - None - stdio handle is inherited.
       - str - opened as a file
       - int -  file number
       - file-like object - fileno() is dupped
       - a PIn or POut derived object
       - a Dev derived object, will create PIn/POut object to link

    If stderr is an instance of DataReader, then stderr is included in
    ProcException on process error.  If the class DataReader is passed
    in as stderr, a DataReader object is created.
    """

    def __init__(self, dag, cmd, stdin=None, stdout=None, stderr=None):
        "setup process. start() must be call to start process"
        self.cmd = tuple(cmd)
        self.dag = dag
        # stdio and argument Dev association
        self.pins = set()
        self.pouts = set()
        self.stdin = self.__stdioAssoc(stdin, "r")
        self.stdout = self.__stdioAssoc(stdout, "w")
        if stderr == DataReader:
            stderr = POut(DataReader())
        self.stderr = self.__stdioAssoc(stderr, "w")
        self.__argsAssoc()
        self.pid = None
        self.statusPipe = None
        self.returncode = None  # exit code, or -signal
        self.exceptInfo = None # (exception, value, traceback)
        self.started = False
        self.finished = False
        self.forced = False    # force termination during ProcDag cleanup

    @staticmethod
    def __devStr(dev):
        "get string used to describe a Dev"
        if isinstance(dev, File):
            return dev.path
        else:
            return str(dev)

    def __str__(self):
        "get simple description of process"
        strs = []
        for a in self.cmd:
            if isinstance(a, PIn):
                strs.append("<(" + self.__devStr(a.dev) + ")")
            elif isinstance(a, PIn):
                strs.append(">(" + self.__devStr(a.dev) + ")")
            else:
                strs.append(_quoteStr(str(a)))
        return " ".join(strs)

    def getPios(self):
        "get set of associated PIn and POut objects"
        return self.pins|self.pouts

    def __stdioAssoc(self, spec, mode):
        """check a stdio spec validity and associate if PInOut or Dev"""
        if (spec is None) or isinstance(spec, int):
            return spec  # passed unchanged
        if callable(getattr(spec, "fileno", None)):
            return spec.fileno()  # is file-like

        # make spec into PInOut object if needed
        if isinstance(spec, str) or isinstance(spec, unicode):
            spec = File(spec)
        if isinstance(spec, Dev):
            if mode == "r":
                spec = PIn(spec)
            else:
                spec = POut(spec)
        elif isinstance(spec, PIn):
            if mode != "r":
                raise Exception("stdout or stderr can not be specified with a PIn object")
        elif isinstance(spec, POut):
            if mode != "w":
                raise Exception("stdin can not be specified with a POut object")
        else:
            raise Exception("invalid stdio specification object type: " + str(type(spec)))

        if isinstance(spec, PIn):
            self.pins.add(spec)
        else:
            self.pouts.add(spec)
        spec.assocByFd(self)
        return spec

    def __argsAssoc(self):
        """call assoc on PInOut arguments"""
        for a in self.cmd:
            if isinstance(a, PIn):
                self.pins.add(a)
                a.assocByPath(self)
            elif isinstance(a, POut):
                self.pouts.add(a)
                a.assocByPath(self)

    def __buildCmd(self):
        """Build command vector, with all object converted to strings.  This
        substitutes references to PInOut objects with strings containing named
        pipe paths."""
        cmd = []
        for a in self.cmd:
            if isinstance(a, PInOut):
                cmd.append(a.getArg())
            else:
                cmd.append(str(a))
        return cmd

    def __stdioSetup(self, spec, stdfd):
        """setup one of the stdio fds."""
        fd = None
        if isinstance(spec, PInOut):
            fd = spec.getFd()
        elif isinstance(spec, str):
            if stdfd == 0:  # stdin?
                fd = os.open(spec, os.O_RDONLY)
            else:
                fd = os.open(spec, os.O_WRONLY|os.O_CREAT|os.O_TRUNC, 0666)
        elif isinstance(spec, int):
            fd = spec
        if (fd is not None) and (fd != stdfd):
            os.dup2(fd, stdfd)
            # Don't close source file here, must delay closing in case stdout/err is same fd

    def __closeFiles(self):
        "clone non-stdio files"
        keepOpen = set([self.statusPipe.wfd]) | trace.getActiveTraceFds()
        for fd in xrange(3, MAXFD+1):
            try:
                if not fd in keepOpen:
                    os.close(fd)
            except: pass

    def __doChildStart(self):
        "guts of start child process"
        self.statusPipe.postForkChild()
        _setPgid(os.getpid(), self.dag.pgid if (self.dag.pgid is not None) else os.getpid())
        cmd = self.__buildCmd()
        self.__stdioSetup(self.stdin, 0)
        self.__stdioSetup(self.stdout, 1)
        self.__stdioSetup(self.stderr, 2)
        self.__closeFiles()
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
        os.execvp(cmd[0], cmd)
            
    def __childStart(self):
        "start in child process"
        try:
            self.__doChildStart()
        except Exception as ex:
            # FIXME: use isinstance(ex, ProcException) causes error in python
            if type(ex) != ProcException:
                ex = ProcException(str(self), cause=ex)
            self.statusPipe.sendExcept(ex)
            
    def __parentStart(self):
        "start in parent process"
        # first process is process leader.
        if self.dag.pgid is None:
            self.dag.pgid = self.pid
        try:
            _setPgid(self.pid, self.dag.pgid)
        except OSError:
            pass # igore error if child has already come and gone
        self.statusPipe.postForkParent()

    def __start(self):
        "do work of starting the process"
        self.statusPipe = _StatusPipe()
        self.started = True  # do first to prevent restarts on error
        self.pid = os.fork()
        if self.pid == 0:
            try:
                self.__childStart()
            finally:
                os.abort() # should never make it here
        else:
            self.__parentStart()

    def _start(self):
        "start the process"
        try:
            self.__start()
        except:
            self.exceptInfo = sys.exc_info()
        if self.exceptInfo is not None:
            self.raiseIfExcept()

    def _execWait(self):
        "receive status, raising the exception if one was send"
        ex = self.statusPipe.recvStatus()
        if isinstance(ex, Exception):
            if not isinstance(ex, ProcException):
                ex = ProcException(str(self), cause=ex)
            raise ex

    def running():
        "determined if this process has been started, but not finished"
        return self.started and not self.finished

    def isRoot(self):
        "a root has no input pipes"
        # don't count pipes to parent
        for pin in self.pins:
            if pin.isPipe() and (pin.dev.pout.proc is not None):
                return False
        return True

    def isLeaf(self):
        "a leaf has no output pipes"
        # don't count pipes to parent
        for pout in self.pouts:
            if pout.isPipe() and (pout.dev.pin.proc is not None):
                return False
        return True

    def raiseIfExcept(self):
        """raise exception if one is saved, otherwise do nothing"""
        if self.exceptInfo is not None:
            raise self.exceptInfo[0], self.exceptInfo[1], self.exceptInfo[2]

    def __handleErrExit(self):
        # get saved stderr, if possible
        stderr = None
        if isinstance(self.stderr, POut) and isinstance(self.stderr.dev, DataReader):
            stderr = self.stderr.dev.get()
        # FIXME: shouldn't save if we killed it
        self.exceptInfo = (ProcException(str(self), self.returncode, stderr), None, None)
        
    def _handleExit(self, waitStat):
        """Handle process exiting, saving status   Call close on all PInOut objects
        to disassociate """
        self.finished = True
        assert(os.WIFEXITED(waitStat) or os.WIFSIGNALED(waitStat))
        self.returncode = os.WEXITSTATUS(waitStat) if os.WIFEXITED(waitStat) else -os.WTERMSIG(waitStat)
        if not ((self.returncode == 0) or (self.returncode == -signal.SIGPIPE)):
            self.__handleErrExit()
        for pin in self.pins:
            pin.close()
        for pout in self.pouts:
            pout.close()

    def _poll(self):
        """Check if all of the processes have completed.  Return True if it
        has, False if it hasn't."""
        if self.finished:
            return True
        w = os.waitpid(self.pid, os.WNOHANG)
        if w[0] != 0:
            self._handleExit(w[1])
        return (w[0] != 0)

    def _forceFinish(self):
        """Forced termination of process.  The forced flag is set, as an
        indication that this was not a primary failure in the pipeline.
        """
        if self.started and not self.finished:
            # check if finished before killing
            if not self._poll():
                self.forced = True
                os.kill(self.pid, signal.SIGKILL)
                w = os.waitpid(self.pid, 0)
                self._handleExit(w[1])

    def failed(self):
        "check if process failed, call after poll() or wait()"
        return (self.exceptInfo is not None)

class _ProcDagDesc(object):
    """Generate a description of a ProcDag for debugging purposes."""
    def __init__(self, dag):
        self.dag = dag
        self.procsSeen = set()
        self.piosSeen = set()   # avoid cycles

    @staticmethod
    def __isPipelinePipe(spec):
        "is spec a PInOut from stdin to stdout"
        if not PInOut.pIsPipe(spec):
            return False
        dev = spec.dev
        if not (PInOut.pHasProc(dev.pin) and PInOut.pHasProc(dev.pout)):
            return False
        return (dev.pout.proc.stdout == dev.pout) and (dev.pin.proc.stdin == dev.pin)

    def __findPipelineStart(self, proc):
        "starting at a proc, walk back stdin->stdout pipeline to process"
        seen = set() # don't hang on cycles
        while self.__isPipelinePipe(proc.stdin) and (not proc in seen):
            seen.add(proc)
            proc = proc.stdin.dev.pout.proc
        return proc

    def __findPipeline(self, proc):
        "find pipeline containing proc, defined by stdin->stdout connections"
        proc = self.__findPipelineStart(proc)
        pline = []
        while not proc in pline:
            pline.append(proc)
            if self.__isPipelinePipe(proc.stdout):
                proc = proc.stdout.dev.pin.proc
            else:
                break
        return pline

    def __partPipelines(self):
        """find linear pipelines and partition into ones whose stdin/out are
        or are not connected to other process args/stderr"""
        notConn = []
        areConn = []
        done = set()
        for proc in self.dag.procs:
            if not proc in done:
                pl = self.__findPipeline(proc)
                done |= set(pl)
                if PInOut.pHasOtherProc(pl[0].stdin) or PInOut.pHasOtherProc(pl[-1].stdout):
                    areConn.append(pl)
                else:
                    notConn.append(pl)
        return (notConn, areConn)

    def __descPInOut(self, pio):
        "generate descriptor of a PInOut"
        if isinstance(pio.dev, Pipe):
            if pio in self.piosSeen:
                return "{CYCLE}"  # hit cycle!
            self.piosSeen.add(pio)
            if isinstance(pio, PIn):
                return self.__descPipeline(pio.dev.pout.proc)
            elif isinstance(pio, POut):
                return self.__descPipeline(pio.dev.pin.proc)
        else:
            return str(pio.dev)

    def __descPInOutArg(self, pio):
        "generate a description of PInOut argument"
        arg = pio.argPrefix if (pio.argPrefix is not None) else ""
        s = self.__descPInOut(pio)
        if isinstance(pio.dev, Pipe):
            if isinstance(pio, PIn):
                sep1 = "<("
            elif isinstance(pio, POut):
                sep1 = ">("
            sep2 = ")"
        else:
            sep1 = sep2 = ""
        return arg + sep1 + s + sep2

    @staticmethod
    def __nonPipeStdioDesc(spec, stdfd, sym):
        "get description to use for a stdio fd when not a pipe"
        if spec is None:
            return ""  # defaulted
        elif isinstance(spec, int):
            if spec != stdfd:
                return " " + sym + "&" + str(spec)
            else:
                return ""  # default, so display nothing
        elif isinstance(spec, PInOut):
            return  " " + sym + str(spec.dev)
        else:
            return  " " + sym + str(spec)

    def __descProc(self, proc):
        """describe a single process in a pipeline, recursively handling args
        and stderr, but not stdin/out in they are pipes"""
        if proc in self.procsSeen:
            # use simplified description if already processed
            return proc.cmd[0] + " ..."
        self.procsSeen.add(proc)
        strs = []
        # command and arguments
        for a in proc.cmd:
            if isinstance(a, PInOut):
                strs.append(self.__descPInOutArg(a))
            else:
                strs.append(_quoteStr(a))
        desc = " " .join(strs)
        # stdin
        if not PInOut.pIsPipe(proc.stdin):
            desc += self.__nonPipeStdioDesc(proc.stdin, 0, "<")
        # stdout
        if not PInOut.pIsPipe(proc.stdout):
            desc += self.__nonPipeStdioDesc(proc.stdout, 1, ">")
        # stderr
        if PInOut.pIsPipe(proc.stderr):
            desc += " 2>(" + self.__descPipeline(proc.stderr.dev.pout.proc) + ")"
        else:
            desc += self.__nonPipeStdioDesc(proc.stderr, 2, "2>")
        return desc

    def __descPipeline(self, proc):
        """depth-first process description generation, starting with the pipeline contain proc"""
        proc = self.__findPipelineStart(proc)
        descs = []
        while True:
            descs.append(self.__descProc(proc))
            if self.__isPipelinePipe(proc.stdout):
                proc = proc.stdout.dev.pin.proc
            else:
                break
        return " | ".join(descs)

    def __str__(self):
        """get a string more or less describing the DAG"""
        # find sub-pipelines not connected as args or stderr and start
        # formatting these
        (notConn, areConn)= self.__partPipelines()
        notConn.sort()  # consistent test results
        areConn.sort()
        descs = []
        for pl in notConn:
            descs.append(self.__descPipeline(pl[0]))
        desc = " ; ".join(descs)

        # handle ones missed, most likely a cycle:
        if len(self.procsSeen) < len(self.dag.procs):
            descs = []
            for p in self.dag.procs:
                if not p in self.procsSeen:
                    descs.append(str(p))
            descs.sort()  # reproducible
            desc += "{CYCLE}: " + " ; ".join(descs)
        return desc
        
class ProcDag(object):
    """Process DAG. Controls creation and management of process graph."""
    def __init__(self):
        self.procs = set()
        self.devs = set()
        self.pgid = None      # process group leader
        self.byPid = dict()   # indexed by pid
        self.started = False  # have procs been started
        self.finished = False # have all procs finished

    def __str__(self):
        """get a string more or less describing the DAG"""
        return str(_ProcDagDesc(self))

    def create(self, cmd, stdin=None, stdout=None, stderr=None):
        "create a new process"
        proc = Proc(self, cmd, stdin, stdout, stderr)
        self.procs.add(proc)
        for pio in proc.getPios():
            self.devs.add(pio.dev)
        return proc

    def getRoots(self):
        "get set of root process, that is, those that don't have input pipes"
        # FIXME: need to think in cliques?
        roots = set()
        for proc in self.procs:
            if proc.isRoot():
                roots.add(proc)
        return roots

    def getLeaves(self):
        "get set of leaf process, that is, those that don't have output pipes"
        # FIXME: needed
        leaves = set()
        for proc in self.procs:
            if proc.isLeaf():
                leaves.add(proc)
        return leaves

    def __dfsValidate(self, proc, seen):
        if proc in seen:
            raise ProcDagException("cycle detected: entering: " + str(proc))
        seen.add(proc)
        for pout in proc.pouts:
            if (pout.dev.pin is not None) and (pout.dev.pin.proc is not None):
                self.__dfsValidate(pout.dev.pin.proc, seen)

    def __validateRoot(self, root):
        seen = set()
        self.__dfsValidate(root, seen)
        return seen

    def __validate(self):
        if len(self.procs) == 0:
            raise ProcDagException("no processes defined")
        roots = self.getRoots()
        if len(roots) == 0:
            raise ProcDagException("cycle detected: no root process")
        seen = set()
        for root in roots:
            seen |= self.__validateRoot(root)
        if seen != self.procs:
            raise ProcDagException("process graph not full connected")

    def __preFork(self):
        for d in self.devs:
            d.preFork()

    def __start(self):
        for p in self.procs:
            p._start()
            self.byPid[p.pid] = p

    def __execBarrier(self):
        for p in self.procs:
            p._execWait()

    def __postExec(self):
        for d in self.devs:
            d.postExecParent()

    def __finish(self):
        "finish up when no errors have occurred"
        self.finished = True
        for d in self.devs:
            d.finish()

    def __cleanupDev(self, dev):
        try:
            dev.finish()
        except Exception as ex:
            # FIXME: make optional, or record, or something
            exi = sys.exc_info()
            stack = "" if exi is None else "".join(traceback.format_list(traceback.extract_tb(exi[2])))+"\n"
            sys.stderr.write("ProcDag dev cleanup exception: " +str(ex)+"\n"+stack)

    def __cleanupProc(self, proc):
        try:
            proc._forceFinish()
        except Exception as ex:
            # FIXME: make optional
            sys.stderr.write("ProcDag proc cleanup exception: " +str(ex)+"\n")
        
    def __cleanup(self):
        """forced cleanup of child processed after failure"""
        self.finished = True
        for d in self.devs:
            self.__cleanupDev(d)
        for p in self.procs:
            self.__cleanupProc(p)

    def start(self):
        """start processes"""
        # FIXME: need non-error here for wait/poll
        self.started = True
        self.__validate()
        # clean up devices and process if there is a failure
        try:
            self.__preFork()
            self.__start()
            self.__execBarrier()
            self.__postExec()
        except:
            self.__cleanup()
            raise

    def raiseIfExcept(self):
        """raise exception if any process has one, otherwise do nothing"""
        for p in self.procs:
            p.raiseIfExcept()

    def poll(self):
        """Check if all of the processes have completed.  Return True if it
        has, False if it hasn't."""
        if not self.started:
            self.start()
        try:
            for p in self.procs:
                if not p._poll():
                    return False
            self.__finish()
        except:
            self.__cleanup()
            raise
        return True

    def __waitOnOne(self):
        "wait on the next process in group to complete, return False if no more"
        try:
            w = os.waitpid(-self.pgid, 0)
        except OSError as ex:
            if ex.errno == errno.ECHILD:
                return False
            raise
        p = self.byPid[w[0]]
        p._handleExit(w[1])
        return True

    def wait(self):
        """Wait for all of the process to complete. Generate an exception if
        any exits non-zero or signals. Starts process if not already
        running."""
        if not self.started:
            self.start()
        try:
            while self.__waitOnOne():
                pass
            self.__finish()
        except:
            self.__cleanup()
            raise
        self.raiseIfExcept()

    def failed(self):
        "check if any process failed, call after poll() or wait()"
        for p in self.procs:
            if p.failed():
                return True
        return False

    def kill(self, sig=signal.SIGTERM):
        "send a signal to the process"
        os.kill(-self.pgid, sig)
        
class Procline(ProcDag):
    """Process pipeline"""
    def __init__(self, cmds, stdin=None, stdout=None, stderr=None):
        """cmds is either a list of arguments for a single process, or a list
        of such lists for a pipeline. If the stdin/out/err arguments are none,
        they are inherited.  Otherwise they can be string file names, a file-like
        object, a file number,  a Dev object, or a PIn or POut object.  Stdin is
        input to the first process, stdout is output to the last process and
        stderr is attached to all processed."""
        ProcDag.__init__(self)
        # FIXME: needed??
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr

        if isinstance(cmds[0], str):
            cmds = [cmds]  # one-process pipeline
        prevPipe = None
        lastCmd = cmds[len(cmds)-1]
        for cmd in cmds:
            prevPipe = self._createProc(cmd, prevPipe, (cmd==lastCmd), stdin, stdout, stderr)
        
    def _createProc(self, cmd, prevPipe, isLastCmd, stdinFirst, stdoutLast, stderr):
        """create one process"""
        if (prevPipe is None):
            stdin = stdinFirst  # first process in pipeline
        else:
            stdin = PIn(prevPipe)
        if (isLastCmd):
            outPipe = None
            stdout = stdoutLast # last process in pipeline
        else:
            outPipe = Pipe()
            stdout = POut(outPipe)
        self.create(cmd, stdin=stdin, stdout=stdout, stderr=stderr)
        return outPipe

class Pipeline(Procline):
    """Object to create and manage a pipeline of processes. It can either run
    an independent set of processes, or a file-like object that either writes
    to or reads from the pipeline.  This provides a simplified interface to
    the ProcDag class, where pipelines are specified as lists.
    """

    # FIXME: change otherEnd stdio, or stdin/stdout, match with mode
    def __init__(self, cmds, mode='r', otherEnd=None):
        """cmds is either a list of arguments for a single process, or
        a list of such lists for a pipeline.  Mode is 'r' for a pipeline
        who's output will be read, or 'w' for a pipeline to that is to
        have data written to it.  If otherEnd is specified, and is a string,
        it is a file to open as stdio file at the other end of the pipeline.
        If it's not a string, it is assumed to be a file object to use for output.
        
        read pipeline ('r'):
          otherEnd --> cmd[0] --> ... --> cmd[n] --> fh
        
        write pipeline ('w')
          fh --> cmd[0] --> ... --> cmd[n] --> otherEnd

        The field fh is the file object used to access the pipeline.
        """
        # FIXME update doc on otherEnd
        if not ((mode == "r") or (mode == "w")):
            raise IOError('invalid mode "' + mode + '"')
        self.mode = mode
        self.closed = False
        self.otherEnd = otherEnd
        self.pio = None

        (otherFh, closeOther) = self._getOtherFh()
        if mode == "r":
            firstIn = otherFh if (otherFh is not None) else 0
            self.pio = PIn(Pipe())
            lastOut = POut(self.pio.dev)
        else:
            lastOut = otherFh if (otherFh is not None) else 1
            self.pio = POut(Pipe())
            firstIn = PIn(self.pio.dev)
        Procline.__init__(self, cmds, stdin=firstIn, stdout=lastOut)
        self.start()
        self.fh = self.pio.getFh()

    def __enter__(self):
        "support for with statement"
        return self

    def __exit__(self, type, value, traceback):
        "support for with statement"
        self.close()

    def _getOtherFh(self):
        """get the other end of the pipeline, return (otherFh, closeOther), with otherFh
        being None if the other end was not opened"""
        if self.otherEnd is None:
            otherFh = None
            closeOther = False
        elif isinstance(self.otherEnd, str):
            if self.mode == 'r':
                otherFh = PIn(File(self.otherEnd))
            else:
                otherFh = POut(File(self.otherEnd))
            closeOther = True
        else:
            otherFh = self.otherEnd
            closeOther = False
        return (otherFh, closeOther)

    def __iter__(self):
        "iter over contents of file"
        return self.fh.__iter__()

    def next(self):
        return self.fh.next()
  
    def flush(self):
        "Flush the internal I/O buffer."
        self.fh.flush()

    def fileno(self):
        "get the integer OS-dependent file handle"
        return self.fh.fileno()
  
    def write(self, str):
        "Write string str to file."
        self.fh.write(str)

    def writeln(self, str):
        "Write string str to file followed by a newline."
        self.fh.write(str)
        self.fh.write("\n")

    def read(self, size=-1):
        return self.fh.read(size)

    def readline(self, size=-1):
        return self.fh.readline(size)

    def readlines(self, size=-1):
        return self.fh.readlines(size)

    def wait(self):
        """wait to for processes to complete, generate an exception if one
        exits no-zero"""
        if self.fh is not None:
            self.pio.close()
            self.fh = None
        Procline.wait(self)

    def poll(self):
        "poll is not allowed for Pipeline objects"
        # don't know what to do about our open pipe, so disallow it
        raise ProcDagException("Pipeline.poll() is not supported")

    def close(self):
        "wait for process to complete, with an error if it exited non-zero"
        if not self.finished:
            self.wait()

__all__ = [ProcException.__name__, PIn.__name__, POut.__name__, Dev.__name__,
           DataReader.__name__, DataWriter.__name__, Pipe.__name__, File.__name__,
           Proc.__name__, ProcDag.__name__, Procline.__name__, Pipeline.__name__]
