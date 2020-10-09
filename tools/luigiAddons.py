"""
Addons for Luigi. Includes decorators multiple inheritance and requirements as well as abstract classes extending
both the Task and Target paradigms.
"""
import sqlite3
import luigi


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
        @luigi.task._task_wraps(task_that_inherits)
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
        @luigi.task._task_wraps(task_that_requires)
        class Wrapped(task_that_requires):
            def requires(self):
                return (self.clone(x) for x in tasks_to_require)

        return Wrapped


class IndexTarget(luigi.Target):
    """
    luigi target that determines if the indices have been built on a hints database.
    """

    def __init__(self, db):
        self.db = db

    def exists(self, timeout=6000):
        con = sqlite3.connect(self.db, timeout=timeout)
        cur = con.cursor()
        r = []
        for idx in ["gidx", "hidx"]:
            query = 'PRAGMA index_info("{}")'.format(idx)
            try:
                v = cur.execute(query).fetchall()
            except sqlite3.OperationalError as exc:
                raise RuntimeError("query failed: {}\nOriginal error message: {}".format(query, exc))
            if len(v) > 0:
                r.append(v)
        return len(r) == 2
