"""
Miscellaneous tools for the pipeline. Some may eventually be refactored into their own modules.
"""
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
