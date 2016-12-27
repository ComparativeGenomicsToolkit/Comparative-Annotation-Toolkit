"""
Miscellaneous tools for the pipeline. Some may eventually be refactored into their own modules.
"""
import itertools
import argparse
import pysam

import procOps
from pipeline import ProcException


class HashableNamespace(argparse.Namespace):
    """
    Adds a __hash__ function to argparse's Namespace. Follows best practices for implementation of __hash__.
    """
    def __hash__(self):
        def xor(x, y):
            return x ^ hash(y)
        val_iter = self.__dict__.itervalues()
        first = hash(val_iter.next())
        return reduce(xor, val_iter, first) ^ hash(tuple(self.__dict__.values()))


class PipelineNamespace(object):
    """
    A Hashable namespace that maintains knowledge of whether a member is significant and thus should be hashed.
    Used to maintain information on the pipeline state but allow users to change insignificant features without forcing
    the pipeline to rerun expensive modules.
    """
    def __init__(self):
        self.significant = {}

    def set(self, name, val, significant=True):
        setattr(self, name, val)
        self.significant[name] = significant

    def __hash__(self):
        def xor(x, y):
            return x ^ hash(y)
        vals = tuple(name for name in self.__dict__ if name != 'significant' and self.significant[name])
        val_iter = iter(vals)
        first = hash(val_iter.next())
        return reduce(xor, val_iter, first) ^ hash(vals)


def convert_gtf_gp(gp_target, gtf_target):
    """converts a GTF to genePred"""
    cmd = ['gtfToGenePred', '-genePredExt', gtf_target.path, '/dev/stdout']
    with gp_target.open('w') as outf:
        procOps.run_proc(cmd, stdout=outf)


def convert_gp_gtf(gtf_target, gp_target, source='CAT'):
    """Converts a genePred to GTF"""
    cmd = ['genePredToGtf', 'file', gp_target.path, '-utr', '-honorCdsStat', '-source={}'.format(source), '/dev/stdout']
    with gtf_target.open('w') as outf:
        procOps.run_proc(cmd, stdout=outf)


def is_exec(program): 
    """checks if a program is in the global path and executable"""
    cmd = ['which', program]
    try:
        return procOps.call_proc_lines(cmd)[0].endswith(program)
    except ProcException:
        return False


def is_bam(path):
    """Checks if a path is a BAMfile"""
    try:
        pysam.Samfile(path)
    except IOError:
        raise RuntimeError('Path {} does not exist'.format(path))
    except ValueError:
        return False
    return True


def pairwise(iterable):
    """s -> (s0, s1), (s2, s3), (s4, s5), ..."""
    a = iter(iterable)
    return itertools.izip(a, a)
