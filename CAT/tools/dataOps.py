"""
Operations on dictionaries and lists.
"""
from collections import namedtuple
import operator
import itertools
import pandas as pd


def combine_dicts(a, b, op=operator.add):
    """
    http://stackoverflow.com/questions/11011756/is-there-any-pythonic-way-to-combine-two-dicts-adding-values-for-keys-that-appe
    """
    return dict(a.items() + b.items() + [(k, op(a[k], b[k])) for k in b.viewkeys() & a.viewkeys()])


def dict_to_named_tuple(d, name):
    """
    Converts a dict to a named tuple, whose name is name.
    """
    return namedtuple(name, d.keys())(**d)


def merge_dicts(list_of_dicts):
    """
    This will merge a list of dicts. Any duplicate keys will end up with the last value seen.
    """
    return reduce(lambda a, d: a.update(d) or a, list_of_dicts, {})


def flatten_list_of_lists(l):
    """
    http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
    """
    return [item for sublist in l for item in sublist]


def flatten_defaultdict_list(d):
    return {k: flatten_list_of_lists(v) for k, v in d.iteritems()}


def grouper(iterable, size):
    """
    http://stackoverflow.com/questions/434287/what-is-the-most-pythonic-way-to-iterate-over-a-list-in-chunks
    """
    it = iter(iterable)
    chunk = tuple(itertools.islice(it, size))
    while chunk:
        yield chunk
        chunk = tuple(itertools.islice(it, size))


def munge_nested_dicts_for_plotting(data_dict, norm=False, sort_column=None):
    """
    Munges nested dictionaries into a pandas DataFrame. If sort_column is not None, will order rows based on values
    in that column. If sort_column is None, will maintain internal dictionary orders. Obviously, if the internal
    dictionaries are not ordered, then this will be meaningless. This is necessary because while
    pd.DataFrame.from_dict() honors an outer OrderedDict, it does not honor nested ones.
    If norm is True, will normalize so that the sum of every column is 1.
    """
    df = pd.DataFrame.from_dict(data_dict)
    if sort_column is not None:
        assert sort_column < df.shape[1], "sort_column larger than number of columns"
        first_column = df[df.columns[0]].copy()
        first_column.sort(ascending=False)
        row_order = first_column.index
    else:
        row_order = data_dict[data_dict.keys()[0]].keys()
    df = df.reindex(row_order)
    df = df.fillna(0)
    if norm is True:
        df = 100 * df.div(df.sum(axis=0), axis=1)
    return [[x[0], x[1].tolist()] for x in df.iteritems()], list(row_order)
