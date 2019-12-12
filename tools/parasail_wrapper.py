"""
Utility functions that add functionality to the parasail pairwise alignment library.

"""
import re
import parasail
from .misc import pairwise
from .psl import PslRow

cigar_re = re.compile('([MIDNSHPX=])')
INS = 'I'
DEL = 'D'
MATCH = '='
MISMATCH = 'X'


def iter_cigar(cigar):
    ref_pos = cigar.beg_ref
    tgt_pos = cigar.beg_query
    for num, op in pairwise(re.split(cigar_re, cigar.decode)):
        num = int(num)
        yield ref_pos, tgt_pos, num, op
        if op == MATCH or op == MISMATCH:
            ref_pos += num
            tgt_pos += num
        elif op == DEL:
            tgt_pos += num
        elif op == INS:
            ref_pos += num
        else:
            assert False


def construct_fa(seq1, seq2, cigar):
    aln1 = []
    aln2 = []
    for ref_pos, tgt_pos, num, op in iter_cigar(cigar):
        if op == MATCH or op == MISMATCH:
            aln1.append(seq1[ref_pos:ref_pos + num])
            aln2.append(seq2[tgt_pos:tgt_pos + num])
        elif op == DEL:
            aln1.append(''.join(['-'] * min(num, len(seq2) - tgt_pos)))
            aln2.append(seq2[tgt_pos:tgt_pos + num])
        elif op == INS:
            aln1.append(seq1[ref_pos:ref_pos + num])
            aln2.append(''.join(['-'] * min(num, len(seq1) - ref_pos)))
        assert len(aln1[-1]) == len(aln2[-1])
    aln1 = ''.join(aln1)
    aln2 = ''.join(aln2)
    assert len(aln1) == len(aln2)
    assert max(len(seq1), len(seq2)) == len(aln1)
    return aln1, aln2


def construct_psl(name1, name2, result):
    block_sizes = []
    q_starts = []
    t_starts = []

    matches = 0
    mismatches = 0
    q_num_insert = 0
    q_base_insert = 0
    t_num_insert = 0
    t_base_insert = 0
    q_pos = result.cigar.beg_query
    q_size = result.len_query
    t_pos = result.cigar.beg_ref
    t_size = result.len_ref

    parsed_cigar = list(pairwise(re.split(cigar_re, result.cigar.decode.decode('utf-8'))))

    for i, (num, op) in enumerate(parsed_cigar):
        num = int(num)
        if op == MATCH or op == MISMATCH:
            block_sizes.append(num)
            q_starts.append(q_pos)
            t_starts.append(t_pos)
            t_pos += num
            q_pos += num
            if op == MATCH:
                matches += num
            else:
                mismatches += num
        elif op == INS:
            q_pos += num
            q_num_insert += 1
            q_base_insert += num
        elif op == DEL:
            # ignore deletions on the ends
            if i == 0 or i == len(parsed_cigar) - 1:
                continue
            t_pos += num
            t_num_insert += 1
            t_base_insert += num
        else:
            assert False
    block_count = len(block_sizes)
    p = PslRow((matches, mismatches, 0, 0, q_num_insert, q_base_insert, t_num_insert, t_base_insert, '+',
                name1, q_size, q_starts[0], q_starts[block_count - 1] + block_sizes[block_count - 1],
                name2, t_size, t_starts[0], t_starts[block_count - 1] + block_sizes[block_count - 1],
                block_count, ','.join(map(str, block_sizes)), ','.join(map(str, q_starts)), ','.join(map(str, t_starts))))
    return p


def aln_proteins(seq1, name1, seq2, name2):
    result = parasail.sg_trace_scan_32(seq1, seq2, 10, 1, parasail.blosum62)
    return construct_psl(name1, name2, result)


def aln_nucleotides(seq1, name1, seq2, name2):
    result = parasail.sg_trace_scan_32(seq1, seq2, 10, 1, parasail.nuc44)
    return construct_psl(name1, name2, result)
