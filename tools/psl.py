"""
Convenience library for working with psl alignment files.

Original Author: Dent Earl
Modified by Ian Fiddes
"""
import math
from collections import Counter, defaultdict

from tools.fileOps import iter_lines
from tools.mathOps import format_ratio
from tools.nameConversions import strip_alignment_numbers

__author__ = 'Ian Fiddes'


class PslRow(object):
    """ Represents a single row in a PSL file.
    http://genome.ucsc.edu/FAQ/FAQformat.html#format2
    """
    __slots__ = ('matches', 'mismatches', 'repmatches', 'n_count', 'q_num_insert', 'q_base_insert', 't_num_insert',
                 't_base_insert', 'strand', 'q_name', 'q_size', 'q_start', 'q_end', 't_name', 't_size', 't_start',
                 't_end', 'block_count', 'block_sizes', 'q_starts', 't_starts')

    def __init__(self, data_tokens):
        assert(len(data_tokens) == 21)
        self.matches = int(data_tokens[0])
        self.mismatches = int(data_tokens[1])
        self.repmatches = int(data_tokens[2])
        self.n_count = int(data_tokens[3])
        self.q_num_insert = int(data_tokens[4])
        self.q_base_insert = int(data_tokens[5])
        self.t_num_insert = int(data_tokens[6])
        self.t_base_insert = int(data_tokens[7])
        self.strand = data_tokens[8]
        self.q_name = data_tokens[9]
        self.q_size = int(data_tokens[10])
        self.q_start = int(data_tokens[11])
        self.q_end = int(data_tokens[12])
        self.t_name = data_tokens[13]
        self.t_size = int(data_tokens[14])
        self.t_start = int(data_tokens[15])
        self.t_end = int(data_tokens[16])
        self.block_count = int(data_tokens[17])
        # lists of ints
        self.block_sizes = [int(x) for x in data_tokens[18].split(',') if x]
        self.q_starts = [int(x) for x in data_tokens[19].split(',') if x]
        self.t_starts = [int(x) for x in data_tokens[20].split(',') if x]

    def target_coordinate_to_query(self, p):
        """ Take position P in target coordinates (positive) and convert it
        to query coordinates (positive).
        """
        if self.strand not in ['+', '-', '++']:
            raise NotImplementedError('PslRow does not support coordinate conversions for strand {}'.format(self.strand))
        if p < self.t_start:
            return None
        if p >= self.t_end:
            return None
        for i, t in enumerate(self.t_starts):
            if p < t:
                continue
            if p >= t + self.block_sizes[i]:
                continue
            # p must be in block
            offset = p - t
            if self.strand == '+' or self.strand == '++':
                return self.q_starts[i] + offset
            else:
                return self.q_size - (self.q_starts[i] + offset) - 1
        return None

    def query_coordinate_to_target(self, p):
        """ Take position P in query coordinates (positive) and convert it
        to target coordinates (positive).
        """
        if self.strand not in ['+', '-', '++']:
            raise NotImplementedError('PslRow does not support coordinate conversions for strand {}'.format(self.strand))
        if p < self.q_start:
            return None
        if p >= self.q_end:
            return None
        if self.strand == '-':
            p = self.q_size - p - 1
        for i, q in enumerate(self.q_starts):
            if p < q:
                continue
            if p >= q + self.block_sizes[i]:
                continue
            # p must be in block
            offset = p - q
            return self.t_starts[i] + offset
        return None

    @property
    def coverage(self):
        return format_ratio(self.matches + self.mismatches + self.repmatches, self.q_size,
                            num_digits=5, resolve_nan=0)

    @property
    def identity(self):
        return format_ratio(self.matches + self.repmatches,
                            self.matches + self.repmatches + self.mismatches + self.q_base_insert,
                            num_digits=5, resolve_nan=0)

    @property
    def target_identity(self):
        return format_ratio(self.matches + self.repmatches,
                            self.matches + self.repmatches + self.mismatches + self.q_base_insert + self.t_base_insert,
                            num_digits=5, resolve_nan=0)

    @property
    def target_coverage(self):
        return format_ratio(self.matches + self.mismatches + self.repmatches, self.t_size,
                            num_digits=5, resolve_nan=0)

    @property
    def percent_n(self):
        return format_ratio(self.n_count, self.q_size, num_digits=5)

    @property
    def badness(self):
        """
        The Jim Kent Badness score, attempts to calculate how bad this alignment is

        https://github.com/ucscGenomeBrowser/kent/blob/fb80e018778062c49021f2c35607868df1054e8e/src/hg/pslCDnaFilter/cDnaAligns.c#L52-L70

        (mismatch + ref_insert + 3 * log(1 + max(ref_size - tgt_size, 0))) / (match + mismatch + repmatch)

        :return: A float between 0 and 1 where 1 is very bad
        """
        b = format_ratio(self.mismatches + self.q_num_insert + 3 * math.log(1 + abs(self.q_size - self.t_size)),
                         self.matches + self.mismatches + self.repmatches,
                         num_digits=5, resolve_nan=1)
        return min(b, 1)

    def psl_string(self):
        """
        Return a list capable of producing a new PslRow object
        """
        return map(str, [self.matches, self.mismatches, self.repmatches, self.n_count, self.q_num_insert,
                         self.q_base_insert, self.t_num_insert, self.t_base_insert, self.strand, self.q_name,
                         self.q_size, self.q_start, self.q_end, self.t_name, self.t_size, self.t_start,
                         self.t_end, self.block_count, ','.join([str(b) for b in self.block_sizes]),
                         ','.join([str(b) for b in self.q_starts]),
                         ','.join([str(b) for b in self.t_starts])])


def psl_iterator(psl_file, make_unique=False):
    """
    Iterates over PSL file generating PslRow objects returning the name and the object itself
    """
    counts = Counter()
    with open(psl_file) as inf:
        for tokens in iter_lines(inf):
            psl = PslRow(tokens)
            if make_unique is True:
                numbered_aln_id = '-'.join([psl.q_name, str(counts[psl.q_name])])
                counts[psl.q_name] += 1
                psl.q_name = numbered_aln_id
            yield psl


def get_alignment_dict(psl_file, make_unique=False):
    """
    Convenience function for creating a dictionary of PslRow objects.
    """
    return {psl.q_name: psl for psl in psl_iterator(psl_file, make_unique)}


def group_alignments_by_qname(psl_iterator, strip=True):
    """Convenience function for grouping PSLs by their name"""
    r = defaultdict(list)
    if strip is True:
        for p in psl_iterator:
            r[strip_alignment_numbers(p.q_name)].append(p)
    else:
        for p in psl_iterator:
            r[p.q_name].append(p)
    return r
