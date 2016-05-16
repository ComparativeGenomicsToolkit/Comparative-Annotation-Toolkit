"""
Abstract classes for use by luigi.
"""
import luigi
from tools.procOps import run_proc


class AbstractAtomicFileTask(luigi.Task):
    """
    Abstract Task for single files.
    """
    def run_cmd(self, cmd):
        """
        Run a external command that will produce the output file for this task to stdout. Capture this to the file.
        """
        # luigi localTargets guarantee atomicity if used as a handle
        with self.output().open('w') as outf:
            run_proc(cmd, stdout=outf)

