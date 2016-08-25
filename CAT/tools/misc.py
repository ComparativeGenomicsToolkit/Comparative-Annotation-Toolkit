"""
Miscellaneous tools for the pipeline. Some may eventually be refactored into their own modules.
"""
import pysam
import itertools
import procOps


def convert_gtf_gp(gp_target, gtf_target):
    """converts the Augustus output GTF to genePred"""
    with gp_target.open('w') as outf:
        cmd = ['gtfToGenePred', '-genePredExt', gtf_target.path, '/dev/stdout']
        procOps.run_proc(cmd, stdout=outf)


def is_exec(program): 
    """checks if a program is in the global path and executable"""
    cmd = ['which', program]
    try:
        return procOps.call_proc_lines(cmd)[0].endswith(program)
    except Exception as error:
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
