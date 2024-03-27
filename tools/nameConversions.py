"""
Perform name conversions on transMap/AugustusTMR transcripts.
"""

import re


def remove_alignment_number(aln_id, aln_re=re.compile("-[0-9]+$")):
    """
    If the name of the transcript ends with -d as in
    ENSMUST00000169901.2-1, return ENSMUST00000169901.2
    :param aln_id: name string
    :param aln_re: compiled regular expression
    :return: string
    """
    return aln_re.split(aln_id)[0]


def remove_augustus_alignment_number(aln_id, aug_re=re.compile("^aug(TM|TMR|CGP|PB)-")):
    """
    removes the alignment numbers prepended by AugustusTM/AugustusTMR
    Format: aug(TM|TMR|CGP|PB)-ENSMUST00000169901.2-1
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
    return True if remove_augustus_alignment_number(aln_id) == aln_id and remove_alignment_number(aln_id) != aln_id else False


def aln_id_is_augustus_tm(aln_id):
    return aln_id.startswith('augTM-')


def aln_id_is_augustus_tmr(aln_id):
    return aln_id.startswith('augTMR-')


def aln_id_is_cgp(aln_id):
    return aln_id.startswith('augCGP-')


def aln_id_is_pb(aln_id):
    return aln_id.startswith('augPB-')


def aln_id_is_exref(aln_id):
    return aln_id.startswith('exRef-')


def aln_id_is_denovo(aln_id):
    return aln_id_is_pb(aln_id) or aln_id_is_cgp(aln_id)


def alignment_type(aln_id):
    """returns what type of alignment this ID is"""
    if aln_id_is_augustus_tmr(aln_id):
        return 'augTMR'
    elif aln_id_is_augustus_tm(aln_id):
        return 'augTM'
    elif aln_id_is_cgp(aln_id):
        return 'augCGP'
    elif aln_id_is_pb(aln_id):
        return 'augPB'
    elif aln_id_is_exref(aln_id):
        return 'exRef'
    elif aln_id_is_transmap(aln_id):
        return 'transMap'
