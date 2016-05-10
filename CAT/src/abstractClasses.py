"""
Abstract classes for use by luigi.
"""
import luigi
import shutil
import argparse
from toil.job import Job
from tools.procOps import run_proc
from tools.fileOps import atomic_install, ensure_file_dir, ensure_dir


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


class AbstractAtomicManyFileTask(luigi.Task):
    """
    Abstract Task for many files. Used if a program outputs multiple files or cannot write to stdout.
    """
    def get_tmp(self):
        return luigi.LocalTarget(is_tmp=True)

    def run_cmd(self, cmd, tmp_files):
        """
        Run a external command that will produce the output file for this task to many files.
        These files will be atomically installed.
        """
        run_proc(cmd)
        for tmp_f, f in zip(*(tmp_files, self.output())):
            f.makedirs()
            if isinstance(tmp_f, luigi.LocalTarget):
                atomic_install(tmp_f.path, f.path)
            elif isinstance(tmp_f, str):
                atomic_install(tmp_f, f.path)
            else:
                raise NotImplementedError
