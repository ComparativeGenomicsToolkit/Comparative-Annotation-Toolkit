"""
Library mathematical operations
"""
import bisect
import math
from operator import itemgetter
from itertools import groupby


def format_ratio(numerator, denominator, num_digits=None, resolve_nan=None):
    """
    Convenience function that converts two numbers to a ratio.
    Handles dividing by zero, as well as transforming values into floats.
    Rounds the number to the number of num_digits, if requested (not None)
    resolve_nan defines what to do when dividing by zero. Default is to return float('nan'), but this can be changed.
    """
    if denominator == 0 or math.isnan(denominator) or math.isnan(numerator):
        if resolve_nan is None:
            return float('nan')
        else:
            return resolve_nan
    r = float(numerator) / float(denominator)
    if num_digits is not None:
        r = round(r, num_digits)
    return r


def find_closest(numeric_list, query_number):
    """
    Given a list of numbers, and a single query number, find the number in the sorted list that is numerically
    closest to the query number. Uses list bisection to do so, and so should be O(log n)
    """
    sorted_numeric_list = sorted(numeric_list)
    pos = bisect.bisect_left(sorted_numeric_list, query_number)
    if pos == 0:
        return sorted_numeric_list[0]
    if pos == len(sorted_numeric_list):
        return sorted_numeric_list[-1]
    before = sorted_numeric_list[pos - 1]
    after = sorted_numeric_list[pos]
    if after - query_number < query_number - before:
        return after
    else:
        return before


def all_disjoint(sets):
    """http://stackoverflow.com/questions/22432814/check-if-a-collection-of-sets-is-pairwise-disjoint"""
    all = set()
    for s in sets:
        for x in s:
            if x in all:
                return False
            all.add(x)
    return True


def find_intervals(data):
    for k, g in groupby(enumerate(data), lambda (i, x): i - x):
        yield map(itemgetter(1), g)
