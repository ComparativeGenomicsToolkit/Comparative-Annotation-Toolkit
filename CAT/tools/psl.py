"""
Convenience library for working with psl alignment files.

Original Author: Dent Earl
Modified by Ian Fiddes
"""
import re
from collections import Counter
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


def psl_iterator(psl_file, make_unique=False):
    """
    Iterates over PSL file generating PslRow objects returning the name and the object itself
    """
    counts = Counter()
    with open(psl_file) as inf:
        for tokens in iter_lines(inf):
            psl = PslRow(tokens)
            if make_unique is True:
                counts[psl.q_name] += 1
                numbered_aln_id = '-'.join([psl.q_name, str(counts[psl.q_name])])
                psl.q_name = numbered_aln_id
            yield psl.q_name, psl


def get_alignment_dict(psl_file, make_unique=False):
    """
    Convenience function for creating a dictionary of PslRow objects.
    """
    if make_unique is False:
        return {aln_id: aln for aln_id, aln in psl_iterator(psl_file, make_unique)}


def remove_alignment_number(aln_id, aln_re=re.compile("-[0-9]+$")):
    """
    If the name of the transcript ends with -d as in
    ENSMUST00000169901.2-1, return ENSMUST00000169901.2
    :param aln_id: name string
    :param aln_re: compiled regular expression
    :return: string
    """
    return aln_re.split(aln_id)[0]


def remove_augustus_alignment_number(aln_id, aug_re=re.compile("^((augI[0-9]+-[0-9]+)|(augI[0-9]+))-")):
    """
    removes the alignment numbers prepended by augustus
    :param aln_id: name string
    :param aug_re: compiled regular expression
    :return: string
    """
    return aug_re.split(aln_id)[-1]


def strip_alignment_numbers(aln_id):
    """
    Convenience function for stripping both Augustus and transMap alignment IDs from a aln_id
    :param aln_id: name string
    :return: string
    """
    return remove_alignment_number(remove_augustus_alignment_number(aln_id))


def aln_id_is_augustus(aln_id):
    """
    Uses remove_augustus_alignment_number to determine if this transcript is an Augustus transcript
    :param aln_id: name string
    :return: boolean
    """
    return True if remove_augustus_alignment_number(aln_id) != aln_id else False


def aln_id_is_transmap(aln_id):
    """
    Uses remove_augustus_alignment_number to determine if this transcript is an Augustus transcript
    :param aln_id: name string
    :return: boolean
    """
    return True if remove_alignment_number(aln_id) != aln_id else False


def fix_ref_q_starts(ref_psl):
    """
    Inverts a negative strand reference psl. Needed for fuzzy intron determination.
    :param ref_psl: PslRow object generated by GenePredToFakePsl
    :return: list
    """
    if ref_psl.strand == '-':
        ref_starts = [ref_psl.q_size - (ref_psl.q_starts[i] + ref_psl.block_sizes[i]) for i in
                      xrange(len(ref_psl.q_starts))]
    else:
        ref_starts = ref_psl.q_starts
    return ref_starts


def is_fuzzy_intron(intron, tm_psl, ref_starts, fuzz_distance=12):
    """
    Determines if a intron is within fuzz distance of its aligned partner.
    :param intron: ChromosomeInterval for this intron
    :param tm_psl: PslRow object for the relationship between tm_tx and ref_tx
    :param ref_starts: list of transcript coordinates that are intron boundaries in the reference transcript
    :param fuzz_distance: max distance allowed to be moved in transcript coordinate space
    :return: boolean
    """
    q_gap_start = tm_psl.target_coordinate_to_query(intron.start - 1)
    q_gap_stop = tm_psl.target_coordinate_to_query(intron.stop)
    fuzzed_start = q_gap_start - fuzz_distance
    fuzzed_stop = q_gap_stop + fuzz_distance
    r = [fuzzed_start <= ref_gap <= fuzzed_stop for ref_gap in ref_starts]
    return True if any(r) else False


def is_original_cds_stop(tm_tx, ref_tx, tm_psl):
    """
    Does this transMap transcript have its original CDS stop?
    :param tm_tx: GenePredTranscript object for transMap transcript
    :param ref_tx: GenePredTranscript object for reference transcript
    :param tm_psl: PslRow object for the relationship between tm_tx and ref_tx
    :return: boolean
    """
    for i in xrange(tm_tx.cds_size - 4, tm_tx.cds_size - 1):
        p = tm_tx.chromosome_coordinate_to_cds(tm_psl.query_coordinate_to_target(ref_tx.cds_coordinate_to_transcript(i)))
        if p is None:
            return False
    return True


def is_original_cds_start(tm_tx, ref_tx, tm_psl):
    """
    Does this transMap transcript have its original CDS start?
    :param tm_tx: GenePredTranscript object for transMap transcript
    :param ref_tx: GenePredTranscript object for reference transcript
    :param tm_psl: PslRow object for the relationship between tm_tx and ref_tx
    :return: boolean
    """
    for i in xrange(3):
        p = tm_tx.chromosome_coordinate_to_cds(tm_psl.query_coordinate_to_target(ref_tx.cds_coordinate_to_transcript(i)))
        if p is None:
            return False
    return True


def is_original_tss(tm_psl, tm_tx, fuzz_distance=50):
    """
    Is this transMap alignment +/- fuzz_distance from the original transcription start point?
    :param tm_psl: PslRow object for the relationship between tm_tx and ref_tx
    :param tm_tx: GenePredTranscript object for transMap transcript
    :param fuzz_distance: integer distance
    :return: boolean
    """
    # convert reference tss to target
    p = tm_psl.query_coordinate_to_target(0)
    if p is None:
        return False
    # convert this to transcript space
    tgt_p = tm_tx.chromosome_coordinate_to_mrna(p)
    return 0 - fuzz_distance <= tgt_p <= fuzz_distance


def is_original_tts(tm_psl, tm_tx, fuzz_distance=50):
    """
    Is this transMap alignment +/- fuzz_distance from the original transcription termination point?
    :param tm_psl: PslRow object for the relationship between tm_tx and ref_tx
    :param tm_tx: GenePredTranscript object for transMap transcript
    :param fuzz_distance: integer distance
    :return: boolean
    """
    # convert reference tss to target
    p = tm_psl.query_coordinate_to_target(len(tm_tx) - 1)
    if p is None:
        return False
    # convert this to transcript space
    tgt_p = tm_tx.chromosome_coordinate_to_mrna(p)
    return len(tm_tx) - fuzz_distance <= tgt_p <= len(tm_tx) + fuzz_distance
