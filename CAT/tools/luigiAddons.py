"""
Addons for Luigi. Includes decorators multiple inheritance and requirements as well as abstract classes extending
both the Task and Target paradigms.
"""
import luigi
import luigi.util
import luigi.contrib.sqla
import procOps


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
            procOps.run_proc(cmd, stdout=outf)


class multiple_inherits(object):
    """
    Task inheritance.

    Usage:

    .. code-block:: python

        class TaskPathA(luigi.Task):
            a = luigi.IntParameter()
            # ...

        class TaskPathB(luigi.Task):
            b = luigi.IntParameter()

        @multiple_inherits(TaskPathA, TaskPathB):
        class MyTask(luigi.Task):
            def requires(self):
               return self.clone_parent()

            def run(self):
               print self.a # this will be defined
               print self.b # this will also be defined
               # ...
    """

    def __init__(self, *tasks_to_inherit):
        super(multiple_inherits, self).__init__()
        self.tasks_to_inherit = tasks_to_inherit

    def __call__(self, task_that_inherits):
        tasks_to_inherit = self.tasks_to_inherit
        for task_to_inherit in tasks_to_inherit:
            for param_name, param_obj in task_to_inherit.get_params():
                if not hasattr(task_that_inherits, param_name):
                    setattr(task_that_inherits, param_name, param_obj)

        # Modify task_that_inherits by subclassing it and adding methods
        @luigi.util.task_wraps(task_that_inherits)
        class Wrapped(task_that_inherits):
            def clone_parent(self, **args):
                task = self.clone(cls=tasks_to_inherit[0])
                for additional_task in tasks_to_inherit[1:]:
                    task = task.clone(cls=additional_task, **args)
                return task
        return Wrapped


class multiple_requires(object):
    """
    Same as @multiple_inherits, but also auto-defines the requires method.
    """

    def __init__(self, *tasks_to_require):
        super(multiple_requires, self).__init__()
        self.inherit_decorator = multiple_inherits(*tasks_to_require)
        self.tasks_to_require = tasks_to_require

    def __call__(self, task_that_requires):
        task_that_requires = self.inherit_decorator(task_that_requires)
        tasks_to_require = self.tasks_to_require

        # Modify task_that_requires by subclassing it and adding methods
        @luigi.util.task_wraps(task_that_requires)
        class Wrapped(task_that_requires):
            def requires(self):
                return (self.clone(x) for x in tasks_to_require)
        return Wrapped
