"""
Miscellaneous tools for the pipeline. Some may eventually be refactored into their own modules.
"""
import re
import itertools
import argparse
import pysam
import pandas as pd
import os

import procOps
from pipeline import ProcException, Procline
from cat.exceptions import ToolMissingException
from distutils.version import StrictVersion


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
    if running_in_container():
        # We assume containerized versions don't need to check if the
        # tools are installed--they definitely are, and calling docker
        # just to run "which" can be surprisingly expensive. But we do
        # check for the presence of Docker or Singularity, since that should take
        # only a few ms.
        binary_mode = os.environ.get('CAT_BINARY_MODE')
        cmd = ['which', binary_mode]
        pl = Procline(cmd, stdin='/dev/null', stdout='/dev/null', stderr='/dev/null')
        try:
            pl.wait()
        except ProcException:
            raise ToolMissingException("{0} not found. Either install {0}, or install CAT's dependencies and use --binary-mode local.".format(binary_mode))
    cmd = ['which', program]
    try:
        return procOps.call_proc_lines(cmd)[0].endswith(program)
    except ProcException:
        return False


def samtools_version():
    """checks the version of samtools installed"""
    try:
        r = procOps.call_proc_lines(['samtools', '--version'])
        if StrictVersion(r[0].split()[1].split('-')[0]) < '1.3':
            raise ToolMissingException('samtools version is not >= 1.3.0')
    except ProcException:
        raise ToolMissingException('samtools is not installed')


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


def sort_gff(input_file, output_file):
    """Sorts a GFF format file by column 1 (chromosome) then column 4(start integer)"""
    cmd = [['sort', '-n', '-k4,4', input_file], ['sort', '-s', '-n', '-k5,5'], ['sort', '-s', '-k1,1']]
    procOps.run_proc(cmd, stdout=output_file)


def parse_gtf_attr_line(attr_line):
    """parse a GTF attributes line"""
    attr_line = [x.split(' ') for x in re.split('; +', attr_line.replace('"', ''))]
    attr_line[-1][-1] = attr_line[-1][-1].rstrip().replace(';', '')
    return dict(attr_line)


def parse_gff_attr_line(attr_line):
    """parse a GFF attributes line"""
    attr_line = [x.split('=') for x in re.split('; *', attr_line.replace('"', ''))]
    attr_line[-1][-1] = attr_line[-1][-1].rstrip().replace(';', '')
    return dict(attr_line)


def slice_df(df, ix):
    """
    Slices a DataFrame by an index, handling the case where the index is missing. Handles the case where a single row
    is returned, thus making it a series.
    """
    try:
        r = df.xs(ix)
        if isinstance(r, pd.core.series.Series):
            return pd.DataFrame([r])
        else:
            return r
    except KeyError:
        return pd.DataFrame()

def running_in_container():
    """
    Is CAT trying to run tools inside containers?
    """
    return os.environ.get("CAT_BINARY_MODE") != "local"
