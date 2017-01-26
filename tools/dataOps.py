"""
Operations on dictionaries and lists.
"""
import operator
import itertools


def combine_dicts(a, b, op=operator.add):
    """
    http://stackoverflow.com/questions/11011756/is-there-any-pythonic-way-to-combine-two-dicts-adding-values-for-keys-that-appe
    """
    return dict(a.items() + b.items() + [(k, op(a[k], b[k])) for k in b.viewkeys() & a.viewkeys()])


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
