"""
Convenience library for working with psl alignment files.

Original Author: Dent Earl
Modified by Ian Fiddes
"""
import re
from tools.mathOps import format_ratio
from tools.fileOps import iter_lines

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

    def hash_key(self):
        """ return a string to use as dict key.
        """
        return '%s_%s_%d_%d' % (self.q_name, self.t_name, self.t_start, self.t_end)

    def target_coordinate_to_query(self, p):
        """ Take position P in target coordinates (positive) and convert it
        to query coordinates (positive).
        """
        if p < self.t_start:
            return None
        if p >= self.t_end:
            return None
        if self.strand not in ['+', '-']:
            raise RuntimeError('Unanticipated strand: %s' % self.strand)
        for i, t in enumerate(self.t_starts):
            if p < t:
                continue
            if p >= t + self.block_sizes[i]:
                continue
            # p must be in block
            offset = p - t
            if self.strand == '+':
                return self.q_starts[i] + offset
            else:
                return self.q_size - (self.q_starts[i] + offset) - 1
        return None

    def query_coordinate_to_target(self, p):
        """ Take position P in query coordinates (positive) and convert it
        to target coordinates (positive).
        """
        if p < self.q_start:
            return None
        if p >= self.q_end:
            return None
        if self.strand not in ['+', '-']:
            raise RuntimeError('Unanticipated strand: %s' % self.strand)
        # this is the easier one to write
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
        return 100 * format_ratio(self.matches + self.mismatches + self.repmatches, self.q_size)

    @property
    def identity(self):
        return 100 * format_ratio(self.matches + self.repmatches,
                                  self.matches + self.repmatches + self.mismatches + self.q_num_insert)

    @property
    def target_coverage(self):
        return 100 * format_ratio(self.matches + self.mismatches + self.repmatches, self.t_size)

    @property
    def percent_n(self):
        return 100 * format_ratio(self.n_count, self.q_size)

    def psl_string(self):
        """ return SELF as a psl formatted line.
        """
        return "\t".join(map(str, [self.matches, self.mismatches, self.repmatches, self.n_count, self.q_num_insert,
                                   self.q_base_insert, self.t_num_insert, self.t_base_insert, self.strand, self.q_name,
                                   self.q_size, self.q_start, self.q_end, self.t_name, self.t_size, self.t_start,
                                   self.t_end, self.block_count, ','.join([str(b) for b in self.block_sizes]),
                                   ','.join([str(b) for b in self.q_starts]),
                                   ','.join([str(b) for b in self.t_starts])]))


def psl_iterator(psl_file):
    """
    Iterates over PSL file generating PslRow objects returning the name and the object itself
    """
    with open(psl_file) as inf:
        for tokens in iter_lines(inf):
            psl = PslRow(tokens)
            yield psl.q_name, psl


def get_alignment_dict(psl_file):
    """
    Convenience function for creating a dictionary of PslRow objects.
    """
    return {aln_id: aln for aln_id, aln in psl_iterator(psl_file)}


def remove_alignment_number(s, aln_re=re.compile("-[0-9]+$")):
    """
    If the name of the transcript ends with -d as in
    ENSMUST00000169901.2-1, return ENSMUST00000169901.2
    """
    return aln_re.split(s)[0]


def remove_augustus_alignment_number(s, aug_re=re.compile("^((augI[0-9]+-[0-9]+)|(augI[0-9]+))-")):
    """
    removes the alignment numbers prepended by augustus
    """
    return aug_re.split(s)[-1]


def strip_alignment_numbers(aln_id):
    """
    Convenience function for stripping both Augustus and transMap alignment IDs from a aln_id
    """
    return remove_alignment_number(remove_augustus_alignment_number(aln_id))


def aln_id_is_augustus(aln_id):
    """
    Uses remove_augustus_alignment_number to determine if this transcript is an Augustus transcript
    """
    return True if remove_augustus_alignment_number(aln_id) != aln_id else False


def aln_id_is_transmap(aln_id):
    """
    Uses remove_augustus_alignment_number to determine if this transcript is an Augustus transcript
    """
    return True if remove_alignment_number(aln_id) != aln_id else False


def invert_q_starts(ref_aln):
    """
    Inverts the strand of a PSL
    """
    if ref_aln.strand == "-":
        ref_starts = [ref_aln.q_size - (ref_aln.q_starts[i] + ref_aln.block_sizes[i]) for i in
                      xrange(len(ref_aln.q_starts))]
    else:
        ref_starts = ref_aln.q_starts
    return ref_starts
