"""
Classes and methods for parsing a pairwise FASTA alignment file. Expects that each pair of FASTA records will be a
pairwise alignment. Used for transcript-space evaluations.
"""
import itertools
import mathOps
import dataOps
import bio


class AlignmmentRecord(object):
    """
    Parses a pairwise FASTA alignment. Has methods to access the records within.
    """
    def __init__(self, ref_id, ref_seq, tgt_id, tgt_seq):
        """
        Calculates alignment metrics from a pairwise FASTA alignment.
        :param ref_id: reference sequence name
        :param ref_seq: string of reference sequence
        :param tgt_id: target sequence name
        :param tgt_seq: string of the target sequence
        :return: an AlignmentMetrics object
        """
        match = 0
        repmatch = 0
        mismatch = 0
        query_insert = 0
        target_insert = 0
        n_count = 0
        for ref_col, tgt_col in itertools.izip(ref_seq, tgt_seq):
            if ref_col == '-':
                query_insert += 1
            elif tgt_col == '-':
                target_insert += 1
            elif ref_col == 'N' or tgt_col == 'N':
                n_count += 1
            elif ref_col == tgt_col and ref_col.islower() or tgt_col.islower():
                repmatch += 1
            elif ref_col == tgt_col:
                match += 1
            else:
                mismatch += 1
        self.match = match
        self.mismatch = mismatch
        self.repmatch = repmatch
        self.n_count = n_count
        self.query_insert = query_insert
        self.target_insert = target_insert
        self.query_size = len(ref_seq)
        self.target_size = len(tgt_seq)
        self.query_name = ref_id
        self.target_name = tgt_id
        self.ref_seq = ref_seq
        self.tgt_seq = tgt_seq

    @property
    def coverage(self):
        """
        Coverage is defined as (matches + repmatches + mismatches + n_count) / (query_size)
        :return: A float between 0-100
        """
        return 100.0 * mathOps.format_ratio(self.match + self.repmatch + self.mismatch + self.n_count,
                                            self.query_size)

    @property
    def identity(self):
        """
        Identity is defined as (matches + repmatches) / (matches + mismatches + repmatches)
        :return: A float between 0-100
        """
        return 100.0 * mathOps.format_ratio(self.match + self.repmatch,
                                            self.match + self.mismatch + self.repmatch)
    
    @property
    def goodness(self):
        """
        Related to the Jim Kent Badness score, attempts to calculate how bad this alignment is by looking at how
        much of it is insertions

        100 * (mismatch + query_insert + target_insert) / query_size

        :return: A float between 0-100
        """
        return 100 - 100.0 * mathOps.format_ratio(self.mismatch + self.query_insert + self.target_insert,
                                                  self.query_size)


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
    :return: iterable of (ref_id, ref_seq, tgt_id, tgt_seq) tuples
    """
    for (ref_id, ref_seq), (tgt_id, tgt_seq) in dataOps.grouper(bio.read_fasta(fasta_path), 2):
        yield AlignmmentRecord(ref_id, ref_seq, tgt_id, tgt_seq)


def get_alignment_record_dict(fasta_path):
    """
    Wrapper for parse_paired_fasta that returns it as a dictionary
    :param fasta_path: Path to transcript FASTA
    :return: dictionary of {(ref_id, tgt_id): AlignmentRecord)}
    """
    r = {}
    for aln_rec in parse_paired_fasta(fasta_path):
        assert (aln_rec.query_name, aln_rec.target_name) not in r
        r[(aln_rec.query_name, aln_rec.target_name)] = aln_rec
    return r
