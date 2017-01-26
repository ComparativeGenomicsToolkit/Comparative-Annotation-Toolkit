import traceback
import sys


class PycbioException(Exception):
    """Base class for exceptions.  This implements exception chaining and
    stores a stack trace.

    To chain an exception
       try:
          ...
       except Exception as ex:
          raise PycbioException("more stuff", ex)
    """
    def __init__(self, msg, cause=None):
        """Constructor."""
        if (cause is not None) and (not isinstance(cause, PycbioException)):
            # store stack trace in other Exception types
            exi = sys.exc_info()
            if exi is not None:
                setattr(cause, "stackTrace", traceback.format_list(traceback.extract_tb(exi[2])))
        Exception.__init__(self, msg)
        self.msg = msg
        self.cause = cause
        self.stackTrace = traceback.format_list(traceback.extract_stack())[0:-1]

    def __str__(self):
        "recursively construct message for chained exception"
        desc = self.msg
        if self.cause is not None:
            desc += ",\n    caused by: " + self.cause.__class__.__name__ + ": " +  str(self.cause)
        return desc

    def format(self):
        "Recursively format chained exceptions into a string with stack trace"
        return PycbioException.formatExcept(self)

    @staticmethod
    def formatExcept(ex, doneStacks=None):
        """Format any type of exception, handling PycbioException objects and
        stackTrace added to standard Exceptions."""
        desc = type(ex).__name__ + ": "
        # don't recurse on PycbioExceptions, as they will include cause in message
        if isinstance(ex, PycbioException):
            desc += ex.msg + "\n"
        else:
            desc +=  str(ex) +  "\n"
        st = getattr(ex, "stackTrace", None)
        if st is not None:
            if doneStacks is None:
                doneStacks = set()
            for s in st:
                if s not in doneStacks:
                    desc += s
                    doneStacks.add(s)
        ca = getattr(ex, "cause", None)
        if ca is not None:
            desc += "caused by: " + PycbioException.formatExcept(ca, doneStacks)
        return desc
