"""
Miscellaneous tools for the pipeline. Some may eventually be refactored into their own modules.
"""
import pysam
import procOps


def convert_gtf_gp(out_target, gtf_target):
    """converts the Augustus output GTF to genePred"""
    with out_target.open('w') as outf:
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
