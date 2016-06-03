"""
Converts transMap genePred entries into Augustus hints
"""
import intervals
import psl


def construct_start_stop_hints(tm_tx, tm_psl, start_stop_radius=5, tss_tts_radius=10):
    """
    Construct start/stop/tss/tts hints
    :param tm_tx: GenePredTranscript object for transMap transcript
    :param tm_psl: PslRow object for the relationship between tm_tx and ref_tx
    :param start_stop_radius: Radius to extend CDS start/stop hints
    :param tss_tts_radius: Radius to extend tts/tss hints
    :return: list of 4 ChromosomeInterval objects
    """
    hints = []
    if tm_tx.cds_start_stat == 'cmpl':
        start = tm_tx.thick_start - start_stop_radius if tm_tx.strand == '+' else tm_tx.thick_stop - 3 - start_stop_radius
        stop = tm_tx.thick_start + 3 + start_stop_radius if tm_tx.strand == '+' else tm_tx.thick_stop + start_stop_radius
        c = intervals.ChromosomeInterval(tm_tx.chromosome, start, stop, tm_tx.strand,
                                         data={'score': 0, 'name': 'start'})
        hints.append(c)
    if tm_tx.cds_end_stat == 'cmpl':
        start = tm_tx.thick_stop - 3 - start_stop_radius if tm_tx.strand == '+' else tm_tx.thick_start - start_stop_radius
        stop = tm_tx.thick_stop + start_stop_radius if tm_tx.strand == '+' else tm_tx.thick_start + 3 + start_stop_radius
        c = intervals.ChromosomeInterval(tm_tx.chromosome, start, stop, tm_tx.strand,
                                         data={'score': 0, 'name': 'stop'})
        hints.append(c)
    if psl.is_original_tss(tm_psl, tm_tx):
        start = tm_tx.start - tss_tts_radius if tm_tx.strand == '+' else tm_tx.stop - tss_tts_radius
        stop = tm_tx.start + tss_tts_radius if tm_tx.strand == '+' else tm_tx.stop + tss_tts_radius
        c = intervals.ChromosomeInterval(tm_tx.chromosome, start, stop, tm_tx.strand,
                                         data={'score': 0, 'name': 'tss'})
        hints.append(c)
    if psl.is_original_tts(tm_psl, tm_tx):
        start = tm_tx.stop - tss_tts_radius if tm_tx.strand == '+' else tm_tx.start - tss_tts_radius
        stop = tm_tx.stop + tss_tts_radius if tm_tx.strand == '+' else tm_tx.start + tss_tts_radius
        c = intervals.ChromosomeInterval(tm_tx.chromosome, start, stop, tm_tx.strand,
                                         data={'score': 0, 'name': 'tts'})
        hints.append(c)
    return hints


def construct_exon_hints(tm_tx, utrend_cutoff=10):
    """
    Construct exon hints. We generate a hint for every exon interval, moving the first and last by utrend_cutoff only
    if they do not cut into the CDS.
    :param tm_tx: GenePredTranscript object for transMap transcript
    :param utrend_cutoff: Amount to remove from both UTR ends
    :return: list of ChromosomeInterval objects
    """
    hints = []
    for e in tm_tx.exon_intervals:
        start = e.start
        stop = e.stop
        if e == tm_tx.exon_intervals[0]:
            start = min(start + utrend_cutoff, tm_tx.thick_start)
        if e == tm_tx.exon_intervals[-1]:
            stop = max(stop - utrend_cutoff, tm_tx.thick_stop)
        hints.append(intervals.ChromosomeInterval(tm_tx.chromosome, start, stop, tm_tx.strand,
                                                  data={'score': 0, 'name': 'exonpart'}))
    return hints


def construct_intron_hints(tm_tx, ref_psl, tm_psl, max_intron_size=300000):
    """
    Construct intron hints. Each intron that is within fuzz distance of an original intron gets a intron hint.
    :param tm_tx: GenePredTranscript object for transMap transcript
    :param ref_psl: PslRow object for the relationship between the source transcript and genome as made by
    GenePredToFakePsl
    :param tm_psl: PslRow object for the relationship between tm_tx and ref_tx
    :param max_intron_size: Maximum intron size to consider
    :return: list of ChromosomeInterval objects
    """
    hints = []
    ref_starts = psl.fix_ref_q_starts(ref_psl)
    for i in tm_tx.intron_intervals:
        if len(i) > max_intron_size:
            continue
        if psl.is_fuzzy_intron(i, tm_psl, ref_starts):
            hints.append(intervals.ChromosomeInterval(tm_tx.chromosome, i.start, i.stop, tm_tx.strand,
                                                      data={'score': 0, 'name': 'intron'}))
    return hints


def convert_exonpart_to_cdspart_intronpart(hints, tm_tx):
    """
    Converts the list of exonparts hints to cdsparts/intronparts.
    :param hints: List of ChromosomeIntervals produced by construct_exon_hints
    :return: List of ChromosomeIntervals that are pure CDS or pure UTR with the correct name
    """
    cds_interval = intervals.ChromosomeInterval(tm_tx.chromosome, tm_tx.thick_start, tm_tx.thick_stop,
                                                tm_tx.strand)
    utr_intervals = [intervals.ChromosomeInterval(tm_tx.chromosome, tm_tx.start, tm_tx.thick_start, tm_tx.strand),
                     intervals.ChromosomeInterval(tm_tx.chromosome, tm_tx.thick_stop, tm_tx.stop, tm_tx.strand)]
    converted_hints = []
    for interval in hints:
        cds_intersection = interval.intersection(cds_interval)
        if cds_intersection is not None:
            cds_intersection.data = {'score': 0, 'name': 'CDSpart'}
            converted_hints.append(cds_intersection)
        for utr_interval in utr_intervals:
            utr_intersection = interval.intersection(utr_interval)
            if utr_intersection is not None:
                utr_intersection.data = {'score': 0, 'name': 'UTRpart'}
                converted_hints.append(utr_intersection)
    return converted_hints


def generate_fuzzy_hints(converted_hints, exon_part_margin=12):
    """
    Each hint has 3 parts, with the internal core boundaries having the highest score:

           ----------------                                   ----------------     Score 1 hints
           exon_part_margin-----------------------------------exon_part_margin     Score 2 hint
    Exon   -------------------------------------------------------------------

    :param converted_hints: input hints from convert_exonpart_to_cdspart_intronpart
    :param exon_part_margin: buffer size on each end
    :return: List of ChromosomeIntervals that are ready to be turned into hints
    """
    fuzzy_hints = []
    for h in converted_hints:
        core_start = h.start + exon_part_margin
        core_stop = h.stop - exon_part_margin
        if core_start > core_stop:  # this exon is smaller than exon_part_margin * 2
            core_start = core_stop = (core_start + core_stop) / 2
        fuzzy_hints.append(intervals.ChromosomeInterval(h.chromosome, core_start, core_stop, h.strand,
                                                        data={'name': h.data['name'], 'score': 2}))
        if h.start < core_start:
            fuzzy_hints.append(intervals.ChromosomeInterval(h.chromosome, h.start, core_start, h.strand,
                                                            data={'name': h.data['name'], 'score': 1}))
        if h.stop > core_stop:
            fuzzy_hints.append(intervals.ChromosomeInterval(h.chromosome, core_stop, h.stop, h.strand,
                                                            data={'name': h.data['name'], 'score': 1}))
    return fuzzy_hints


def tm_to_hints(tm_tx, ref_psl, tm_psl):
    """
    Converts a genePred transcript to hints parseable by Augustus.

    Note for anyone reading this code: GFF coordinates are 1-based, genePred coordinates are 0-based
    :param tm_tx: GenePredTranscript object for transMap transcript
    :param ref_psl: PslRow object for the relationship between the source transcript and genome as made by
    GenePredToFakePsl
    :param tm_psl: PslRow object for the relationship between tm_tx and ref_tx
    :return: GFF formatted string.
    """

    def interval_to_gff(interval):
        """convert ChromosomeInterval to gff3. expects data to be a dict with a score and a hint name"""
        attributes = ';'.join(['grp={}'.format(tm_tx.name), 'src=T', 'pri=4'])
        return '\t'.join(map(str, [interval.chromosome, 't2h', interval.data['name'], interval.start + 1,
                                   interval.stop + 1, interval.data['score'], interval.strand, '.', attributes]))

    exon_hints = construct_exon_hints(tm_tx)
    converted_hints = convert_exonpart_to_cdspart_intronpart(exon_hints, tm_tx)
    hints = generate_fuzzy_hints(converted_hints)
    hints.extend(construct_intron_hints(tm_tx, ref_psl, tm_psl))
    hints.extend(construct_start_stop_hints(tm_tx, tm_psl))
    gff_hints = [interval_to_gff(i) for i in hints]
    return '\n'.join(gff_hints)
