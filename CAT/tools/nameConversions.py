"""
Perform name conversions on transMap/AugustusTMR transcripts.
"""

import re


def remove_alignment_number(s, aln_re=re.compile("-[0-9]+$")):
    """
    If the name of the transcript ends with -d as in
    ENSMUST00000169901.2-1, return ENSMUST00000169901.2
    """
    return aln_re.split(s)[0]


def remove_augustus_alignment_number(s, aug_re=re.compile("^aug-I[0-9]+-")):
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
