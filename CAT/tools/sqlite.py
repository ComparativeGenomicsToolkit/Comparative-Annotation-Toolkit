"""
Tools to wrap around the native sqlite package. Necessary due to bugs in how pandas interacts with sqlalchemy.
"""
import sqlite3 as sql

__author__ = "Ian Fiddes"


class ExclusiveSqlConnection(object):
    """Context manager for an exclusive SQL connection"""
    def __init__(self, path, timeout=6000):
        self.path = path
        self.timeout = timeout

    def __enter__(self):
        self.con = sql.connect(self.path, timeout=self.timeout, isolation_level="EXCLUSIVE")
        try:
            self.con.execute("BEGIN EXCLUSIVE")
        except sql.OperationalError:
            raise RuntimeError("Database still locked after {} seconds.".format(self.timeout))
        return self.con

    def __exit__(self, exception_type, exception_val, trace):
        self.con.commit()
        self.con.close()


def attach_database(con, path, name):
    """
    Attaches another database found at path to the name given in the given connection.
    """
    con.execute("ATTACH DATABASE '{}' AS {}".format(path, name))


def open_database(path, timeout=6000):
    """opens a database, returning the connection and cursor objects."""
    con = sql.connect(path, timeout=timeout)
    cur = con.cursor()
    return con, cur

