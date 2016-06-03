"""
Represent continuous genomic coordinates. Allows for coordinate arithmetic.
"""
from bio import reverse_complement

__author__ = 'Ian Fiddes'


class ChromosomeInterval(object):
    """
    Represents a continuous genomic interval.
    interval arithmetic adapted from http://code.activestate.com/recipes/576816-interval/
    """
    __slots__ = ('chromosome', 'start', 'stop', 'strand', 'data')

    def __init__(self, chromosome, start, stop, strand, data=None):
        self.chromosome = str(chromosome)
        self.start = int(start)    # 0 based
        self.stop = int(stop)      # exclusive
        assert self.start <= self.stop
        self.strand = strand       # + or -
        self.data = data

    def __len__(self):
        return abs(self.stop - self.start)

    def __hash__(self):
        return (hash(self.chromosome) ^ hash(self.start) ^ hash(self.stop) ^ hash(self.strand) ^
                hash((self.chromosome, self.start, self.stop, self.strand)))

    def __eq__(self, other):
        return (isinstance(other, type(self)) and
                (self.chromosome, self.start, self.stop, self.strand) ==
                (other.chromosome, other.start, other.stop, other.strand))

    def __ne__(self, other):
        return not self == other

    def __gt__(self, other):
        return (isinstance(other, type(self)) and self.chromosome == other.chromosome and
                (self.start, self.stop) > (other.start, other.stop))

    def __ge__(self, other):
        return (isinstance(other, type(self)) and self.chromosome == other.chromosome and
                (self.start, self.stop) >= (other.start, other.stop))

    def __lt__(self, other):
        return (isinstance(other, type(self)) and self.chromosome == other.chromosome and
                (self.start, self.stop) < (other.start, other.stop))

    def __le__(self, other):
        return (isinstance(other, type(self)) and self.chromosome == other.chromosome and
                (self.start, self.stop) <= (other.start, other.stop))

    def __contains__(self, other):
        return self.start <= other < self.stop

    def __add__(self, other):
        if self.strand != other.strand or self.chromosome != other.chromosome:
            return None
        return ChromosomeInterval(self.chromosome, self.start + other.start, self.stop + other.stop, self.strand)

    def __sub__(self, other):
        if self.strand != other.strand or self.chromosome != other.chromosome:
            return None
        return ChromosomeInterval(self.chromosome, self.start - other.start, self.stop - other.stop, self.strand)

    def __repr__(self):
        if self.data is None:
            return "ChromosomeInterval('{}', {}, {}, '{}')".format(self.chromosome, self.start, self.stop, self.strand)
        else:
            return "ChromosomeInterval('{}', {}, {}, '{}', '{}')".format(self.chromosome, self.start, self.stop,
                                                                         self.strand, self.data)

    @property
    def is_null(self):
        if len(self) == 0:
            return True
        return False

    def intersection(self, other):
        """
        Does this interval intersect with another interval?
        :param other: Another ChromosomeInterval object.
        :return: A ChromosomeInterval object representing the intersection, if it exists, otherwise None.
        """
        if self.strand != other.strand or self.chromosome != other.chromosome:
            return None
        if self > other:
            other, self = self, other
        if self.stop <= other.start:
            return None
        if self == other:
            return self
        if self.stop <= other.stop:
            return ChromosomeInterval(self.chromosome, other.start, self.stop, self.strand)
        else:
            return ChromosomeInterval(self.chromosome, other.start, other.stop, self.strand)

    def complement(self, size):
        """
        returns two new ChromosomeIntervals representing the complement of this interval.
        Requires a input chromosome size.
        :param size: Integer genomic size.
        :return: Two ChromosomeInterval objects representing the complement of this interval and the size.
        """
        assert 0 <= len(self) < size
        return [ChromosomeInterval(self.chromosome, 0, self.start, self.strand),
                ChromosomeInterval(self.chromosome, self.stop, size, self.strand)]

    def union(self, other):
        """
        Returns the union of this ChromosomeInterval and another, if it exists.
        :param other: Another ChromosomeInterval object.
        :return: one ChromosomeInterval if a union exists, otherwise None
        """
        if self == other:
            return self
        if self.intersection(other) is not None:
            return self.hull(other)
        return None

    def hull(self, other):
        """
        Returns a new ChromosomeInterval representing the merged interval of these two, regardless of overlap
        I.E. if one interval is [1, 10) and another is [20, 30) this will return [1, 30)
        :param other: Another ChromosomeInterval object.
        :return: one ChromosomeInterval if a hull exists, otherwise None
        """
        if self.chromosome != other.chromosome:
            return None
        if self > other:
            other, self = self, other
        if self.subset(other):
            return ChromosomeInterval(self.chromosome, other.start, other.stop, self.strand)
        elif other.subset(self):
            return ChromosomeInterval(self.chromosome, self.start, self.stop, self.strand)
        return ChromosomeInterval(self.chromosome, self.start, other.stop, self.strand)

    def overlap(self, other, stranded=False):
        """
        Boolean function - does this ChromosomeInterval overlap another?
        :param other: Another ChromosomeInterval object.
        :param stranded: Consider overlaps only on the same strand?
        :return: True if self overlaps other, otherwise False
        """
        if self.chromosome != other.chromosome:
            return False
        if stranded is True and self.strand != other.strand:
            return False
        if self > other:
            other, self = self, other
        return self.stop > other.start

    def subset(self, other, stranded=False):
        """
        Boolean function - is this ChromosomeInterval a subset of another?
        :param other: Another ChromosomeInterval object.
        :param stranded: Consider overlaps only on the same strand?
        :return: True if self overlaps other, otherwise False
        """
        if self.chromosome != other.chromosome:
            return False
        if stranded is True and self.strand != other.strand:
            return False
        return self.start >= other.start and self.stop <= other.stop

    def proper_subset(self, other, stranded=False):
        """
        Boolean function - is this ChromosomeInterval a proper subset of another?
        :param other: Another ChromosomeInterval object.
        :param stranded: Consider overlaps only on the same strand?
        :return: True if self overlaps other, otherwise False
        """
        if self.chromosome != other.chromosome:
            return False
        if stranded is True and self.strand != other.strand:
            return False
        return self.start > other.start and self.stop < other.stop

    def separation(self, other):
        """
        How far in Euclidean distance is this ChromosomeInterval from another?
        :param other: Another ChromosomeInterval object.
        :return: Integer distance
        """
        if self.chromosome != other.chromosome:
            return None
        if self > other:
            other, self = self, other
        if self.stop > other.start:
            return 0
        else:
            return other.start - self.stop

    def symmetric_separation(self, other):
        """
        How far is this ChromosomeInterval from another? Returns two values, for start and stop distances.
        :param other: Another ChromosomeInterval object.
        :return: Integer tuple (left distance, right distance)
        """
        if self.chromosome != other.chromosome:
            return None
        if self > other:
            other, self = self, other
        return other.start - self.start, other.stop - self.stop

    def get_sequence(self, seq_dict, stranded=True):
        """
        Returns the sequence for this ChromosomeInterval. If stranded is True, reverse complements as necessary.
        :param seq_dict: Dictionary-like object with DNA sequences.
        :param stranded: Should we reverse complement negative strand sequences?
        :return: A sequence string.
        """
        if stranded is False or self.strand is True:
            return seq_dict[self.chromosome][self.start:self.stop]
        if self.strand is False:
            return reverse_complement(seq_dict[self.chromosome][self.start:self.stop])


def gap_merge_intervals(intervals, gap):
    """
    Merge gaps between a iterable of ChromosomeIntervals.
    :param intervals: Iterable of ChromosomeIntervals
    :param gap: integer value of gap size to merge.
    :return: List of new ChromosomeIntervals.
    """
    new_intervals = []
    for interval in sorted(intervals):
        if not new_intervals:
            new_intervals.append(ChromosomeInterval(**vars(interval)))
        elif interval.separation(new_intervals[-1]) <= gap:
            new_intervals[-1] = new_intervals[-1].hull(interval)
        else:
            new_intervals.append(ChromosomeInterval(**vars(interval)))
    return new_intervals


def interval_not_intersect_intervals(intervals, interval):
    """
    Determines if one intervals does not overlap an iterable of intervals
    :param intervals: iterable of ChromosomeIntervals
    :param interval: one ChromosomeInterval
    :return: boolean
    """
    intersections = []
    for target_interval in intervals:
        intersections.append(interval.intersection(target_interval))
    return intersections.count(None) == len(intersections)


def interval_not_within_wiggle_room_intervals(intervals, interval, wiggle_room=0):
    """
    Same thing as interval_not_intersect_intervals but looks within wiggle_room bases to count a valid lack of
    intersection. Wiggle room can exist on either side.
    :param intervals: iterable of ChromosomeIntervals
    :param interval: one ChromosomeInterval
    :return: boolean
    """
    try:
        separation = [sum(interval.symmetric_separation(target_interval)) for target_interval in intervals]
    except TypeError:
        return False
    return not any(x <= 2 * wiggle_room for x in separation)  # we allow wiggle on both sides
