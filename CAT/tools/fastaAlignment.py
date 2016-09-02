"""
Classes and methods for parsing a pairwise FASTA alignment file. Expects that each pair of FASTA records will be a
pairwise alignment. Used for transcript-space evaluations.
"""
import collections
import itertools
import mathOps
import dataOps
import bio
import math


class AlignmmentRecord(object):
    """
    Parses a pairwise FASTA alignment. Has methods to access the records within.
    """
    def __init__(self, ref_id, ref_aln, tgt_id, tgt_aln):
        """
        Calculates alignment metrics from a pairwise FASTA alignment.
        :param ref_id: reference sequence name
        :param ref_aln: string of reference sequence
        :param tgt_id: tgt sequence name
        :param tgt_aln: string of the tgt sequence
        :return: an AlignmentMetrics object
        """
        self.ref_name = ref_id
        self.tgt_name = tgt_id
        self.ref_aln = ref_aln
        self.tgt_aln = tgt_aln
        self.ref_seq = ref_aln.replace('-', '')
        self.tgt_seq = tgt_aln.replace('-', '')
        self.ref_size = len(self.ref_seq)
        self.tgt_size = len(self.tgt_seq)
        match, repmatch, mismatch, ref_insert, tgt_insert, n_count = self._calculate_metrics()
        self.match = match
        self.mismatch = mismatch
        self.repmatch = repmatch
        self.n_count = n_count
        self.ref_insert = ref_insert
        self.tgt_insert = tgt_insert
        self.ref_pos_map, self.tgt_pos_map = self._generate_position_map()
        self.ref_inverse_pos_map, self.tgt_inverse_pos_map = self._generate_inverse_position_map()


    @property
    def coverage(self):
        """
        Coverage is defined as (matches + repmatches + mismatches + n_count) / (ref_size)
        :return: A float between 0-100
        """
        return 100.0 * mathOps.format_ratio(self.match + self.repmatch + self.mismatch + self.n_count,
                                            self.ref_size,
                                            num_digits=5)

    @property
    def identity(self):
        """
        Identity is defined as (matches + repmatches) / (matches + mismatches + repmatches)
        :return: A float between 0-100
        """
        return 100.0 * mathOps.format_ratio(self.match + self.repmatch,
                                            self.match + self.mismatch + self.repmatch,
                                            num_digits=5)
    
    @property
    def badness(self):
        """
        The Jim Kent Badness score, attempts to calculate how bad this alignment is

        https://github.com/ucscGenomeBrowser/kent/blob/fb80e018778062c49021f2c35607868df1054e8e/src/hg/pslCDnaFilter/cDnaAligns.c#L52-L70

        (mismatch + ref_insert + 3 * log(1 + max(ref_size - tgt_size, 0))) / (match + mismatch + repmatch)

        :return: A float
        """
        num = self.mismatch + self.ref_insert + 3 * math.log(1 + max(self.ref_size - self.tgt_size, 0))
        denom = self.match + self.mismatch + self.repmatch
        return mathOps.format_ratio(num, denom, num_digits=5)
    
    @property
    def percent_n(self):
        """
        Calculates the percent N of the tgt sequence
        :return: A float between 0-100
        """
        n_count = self.tgt_aln.count('N')
        return 100.0 * mathOps.format_ratio(n_count, self.tgt_size, num_digits=5)

    def _calculate_metrics(self):
        """
        Calculate alignment metrics
        :return: integer values representing match, repmatch, mismatch, ref_insert, tgt_insert, n_count
        """
        match = 0
        repmatch = 0
        mismatch = 0
        ref_insert = 0
        tgt_insert = 0
        n_count = 0
        for ref_col, tgt_col in itertools.izip(self.ref_aln, self.tgt_aln):
            if ref_col == '-':
                ref_insert += 1
            elif tgt_col == '-':
                tgt_insert += 1
            elif ref_col == 'N' or tgt_col == 'N':
                n_count += 1
            elif ref_col == tgt_col and ref_col.islower() or tgt_col.islower():
                repmatch += 1
            elif ref_col == tgt_col:
                match += 1
            else:
                mismatch += 1
        return match, repmatch, mismatch, ref_insert, tgt_insert, n_count

    def _generate_position_map(self):
        """
        Convert the gapped alignment into a map of positions that map the alignment position to the input cDNA position

        :return: Dictionary mapping alignment positions to sequence positions for both sequences in an AlignmentRecord
        """
        pos_map = collections.defaultdict(dict)
        names = [self.ref_name, self.tgt_name]
        alns = [self.ref_aln, self.tgt_aln]
        tgt_is = {n: 0 for n in [self.ref_name, self.tgt_name]}
        for aln_pos, chars in enumerate(itertools.izip(*alns)):
            for name, tgt_i in tgt_is.iteritems():
                pos_map[name][aln_pos] = tgt_i
            for name, c in itertools.izip(*[names, chars]):
                if c != '-':
                    tgt_is[name] += 1
        return pos_map[self.ref_name], pos_map[self.tgt_name]
    
    def _generate_inverse_position_map(self):
        """
        Calculate the positions within the gapped alignment that correspond to cDNA positions
        :return: Dictionary mapping sequence positions to alignment positions for both sequences in an AlignmentRecord
        """
        ref = {y: x for x, y in self.ref_pos_map.iteritems()}
        tgt = {y: x for x, y in self.tgt_pos_map.iteritems()}
        return ref, tgt


def parse_paired_fasta(fasta_path):
    """
    Parse the output of the transcript alignment script.
    This file is of the format:

    >ref_tx
    ATGC
    >tgt_tx
    ATGG
    <newline>
    :param fasta_path: Path to transcript FASTA
    :return: iterable of (ref_id, ref_aln, tgt_id, tgt_aln) tuples
    """
    for (ref_id, ref_aln), (tgt_id, tgt_aln) in dataOps.grouper(bio.read_fasta(fasta_path), 2):
        yield AlignmmentRecord(ref_id, ref_aln, tgt_id, tgt_aln)


def get_alignment_record_dict(fasta_path):
    """
    Wrapper for parse_paired_fasta that returns it as a dictionary
    :param fasta_path: Path to transcript FASTA
    :return: dictionary of {(ref_id, tgt_id): AlignmentRecord)}
    """
    r = {}
    for aln_rec in parse_paired_fasta(fasta_path):
        assert (aln_rec.ref_name, aln_rec.tgt_name) not in r
        r[(aln_rec.ref_name, aln_rec.tgt_name)] = aln_rec
    return r
