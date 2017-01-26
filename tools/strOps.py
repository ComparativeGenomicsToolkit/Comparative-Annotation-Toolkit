# Copyright 2006-2012 Mark Diekhans
"operations on strings"
import re

# matches one or more whitespaces
spaceRe = re.compile("[ \t\n\v\f\r]+")

def hasSpaces(s):
    "test if there are any whitespace characters in a string"
    return spaceRe.search(s) is not None

def splitAtSpaces(s):
    "split a string at one or more contiguous whitespaces"
    return spaceRe.split(s)

def dup(n, s):
    "make a string with n copies of s"
    l = []
    for i in xrange(n):
        l.append(s)
    return "".join(l)

def emptyOrNone(s):
    "is a string empty of None"
    return (s is None) or (len(s) == 0)

def emptyForNone(s):
    "return an empty string if s is None, else s"
    return "" if s is None else s

def noneForEmpty(s):
    "return non if s is a empty string, else s"
    return None if s == "" else s

__all__ = (hasSpaces.__name__, splitAtSpaces.__name__, dup.__name__, emptyForNone.__name__, noneForEmpty.__name__)

