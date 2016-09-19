"""
Represent either BED12 or genePred transcripts as objects. Allows for conversion of coordinates between
chromosome, mRNA and CDS coordinate spaces. Can slice objects into subsets.
"""
from itertools import izip

from bio import reverse_complement, translate_sequence
from fileOps import iter_lines
from intervals import ChromosomeInterval

__author__ = "Ian Fiddes"


class Transcript(object):
    """
    Represent a transcript record from a bed file.
    """
    __slots__ = ('name', 'strand', 'score', 'thick_start', 'rgb', 'thick_stop', 'start', 'stop', 'intron_intervals',
                 'exon_intervals', 'exons', 'block_sizes', 'block_starts', 'block_count', 'chromosome',
                 'interval', 'coding_interval')

    def __init__(self, bed_tokens):
        self.chromosome = bed_tokens[0]
        self.start = int(bed_tokens[1])
        self.stop = int(bed_tokens[2])
        self.name = bed_tokens[3]
        self.score = int(bed_tokens[4])
        self.strand = bed_tokens[5]
        self.thick_start = int(bed_tokens[6])
        self.thick_stop = int(bed_tokens[7])
        self.rgb = bed_tokens[8]
        self.block_count = int(bed_tokens[9])
        self.block_sizes = [int(x) for x in bed_tokens[10].split(",") if x != ""]
        self.block_starts = [int(x) for x in bed_tokens[11].split(",") if x != ""]
        self.exon_intervals = self._get_exon_intervals()
        self.intron_intervals = self._get_intron_intervals()
        self.interval = self._get_interval()
        self.coding_interval = self._get_coding_interval()

    def __len__(self):
        return sum(len(x) for x in self.exon_intervals)

    def __hash__(self):
        return (hash(self.chromosome) ^ hash(self.start) ^ hash(self.stop) ^ hash(self.strand) ^
                hash((self.chromosome, self.start, self.stop, self.strand)))

    @property
    def cds_size(self):
        """calculates the number of coding bases"""
        l = 0
        for e in self.exon_intervals:
            if self.thick_start < e.start and e.stop < self.thick_stop:
                # squarely in the CDS
                l += e.stop - e.start
            elif e.start <= self.thick_start < e.stop < self.thick_stop:
                # thickStart marks the start of the CDS
                l += e.stop - self.thick_start
            elif e.start <= self.thick_start and self.thick_stop <= e.stop:
                # thickStart and thickStop mark the whole CDS
                l += self.thick_stop - self.thick_start
            elif self.thick_start < e.start < self.thick_stop <= e.stop:
                # thickStop marks the end of the CDS
                l += self.thick_stop - e.start
        return l

    @property
    def num_coding_introns(self):
        """how many coding introns does this transcript have?"""
        return len([i for i in self.intron_intervals if i.subset(self.coding_interval)])

    @property
    def num_coding_exons(self):
        """how many coding exons does this transcript have?"""
        return len([i for i in self.exon_intervals if i.overlap(self.coding_interval)])

    def _get_interval(self):
        """
        Returns a ChromosomeInterval object representing the full span of this transcript.
        """
        return ChromosomeInterval(self.chromosome, self.start, self.stop, self.strand)

    def _get_coding_interval(self):
        """
        Returns a ChromosomeInterval object representing the coding span of this transcript.
        """
        return ChromosomeInterval(self.chromosome, self.thick_start, self.thick_stop, self.strand)

    def _get_exon_intervals(self):
        """
        Builds a list of ChromosomeInterval objects representing the exons of this transcript.
        :return: List of ChromosomeIntervals
        """
        exon_intervals = []
        for block_size, block_start in izip(*(self.block_sizes, self.block_starts)):
            start = self.start + block_start
            stop = self.start + block_start + block_size
            exon_intervals.append(ChromosomeInterval(self.chromosome, start, stop, self.strand))
        return exon_intervals

    def _get_intron_intervals(self):
        """
        Builds a list of ChromosomeInterval objects representing the introns of this transcript.
        :return: List of ChromosomeIntervals
        """
        intron_intervals = []
        for i in xrange(1, len(self.block_starts)):
            stop = self.start + self.block_starts[i]
            start = self.start + self.block_starts[i - 1] + self.block_sizes[i - 1]
            intron_intervals.append(ChromosomeInterval(self.chromosome, start, stop, self.strand))
        return intron_intervals

    def get_bed(self, rgb=None, name=None, new_start=None, new_stop=None):
        """
        Returns BED tokens for this object. Can be sliced into sub regions.
        :param rgb: Set this to modify the RGB field.
        :param name: Set this to modify the name field.
        :param new_start: Set this (in chromosome coordinates) to move the start.
        :param new_stop: Set this (in chromosome coordinates) to move the stop.
        :return: List of values representing a BED entry.
        """
        if new_start is not None and new_stop is not None:
            assert new_start <= new_stop
        if new_start is not None:
            assert new_start >= self.start
        if new_stop is not None:
            assert new_stop <= self.stop
        rgb = self.rgb if rgb is None else rgb
        name = self.name if name is None else name
        if new_start is None and new_stop is None:
            block_starts = ",".join(map(str, self.block_starts))
            block_sizes = ",".join(map(str, self.block_sizes))
            return map(str, [self.chromosome, self.start, self.stop, name, self.score, self.strand,
                             self.thick_start, self.thick_stop, rgb, self.block_count, block_sizes, block_starts])
        elif new_start == new_stop:
            return map(str, [self.chromosome, new_start, new_stop, name, self.score, self.strand,
                             new_start, new_stop, rgb, 1, 0, 0])

        def _move_start(exon_intervals, block_count, block_starts, block_sizes, start, new_start):
            to_remove = len([x for x in exon_intervals if x.start <= new_start and x.stop <= new_start])
            assert to_remove < len(exon_intervals)
            if to_remove > 0:
                block_count -= to_remove
                block_sizes = block_sizes[to_remove:]
                start += block_starts[to_remove]
                new_block_starts = [0]
                for i in xrange(to_remove, len(block_starts) - 1):
                    new_block_starts.append(block_starts[i + 1] - block_starts[i] + new_block_starts[-1])
                block_starts = new_block_starts
            if new_start > start:
                block_sizes[0] += start - new_start
                block_starts[1:] = [x + start - new_start for x in block_starts[1:]]
                start = new_start
            return start, block_count, block_starts, block_sizes

        def _move_stop(exon_intervals, block_count, block_starts, block_sizes, stop, start, new_stop):
            to_remove = len([x for x in exon_intervals if x.stop >= new_stop and x.start >= new_stop])
            assert to_remove < len(exon_intervals)
            if to_remove > 0:
                block_count -= to_remove
                block_sizes = block_sizes[:-to_remove]
                block_starts = block_starts[:-to_remove]
                assert len(block_sizes) == len(block_starts)
                if len(block_sizes) == 0:
                    block_sizes = block_starts = [0]
                    block_count = 1
                stop = start + block_sizes[-1] + block_starts[-1]
            if start + block_starts[-1] < new_stop < stop:
                block_sizes[-1] = new_stop - start - block_starts[-1]
                stop = new_stop
            return stop, block_count, block_starts, block_sizes

        start = int(self.start)
        stop = int(self.stop)
        thick_start = int(self.thick_start)
        thick_stop = int(self.thick_stop)
        exon_intervals = self.exon_intervals
        block_count = int(self.block_count)
        block_starts = map(int, self.block_starts)
        block_sizes = map(int, self.block_sizes)

        if new_start is not None and new_start > self.start:
            start, block_count, block_starts, block_sizes = _move_start(exon_intervals, block_count,
                                                                        block_starts, block_sizes, start,
                                                                        new_start)
        if new_stop is not None and new_stop < self.stop:
            stop, block_count, block_starts, block_sizes = _move_stop(exon_intervals, block_count,
                                                                      block_starts, block_sizes, stop,
                                                                      start, new_stop)
        if start > thick_start:
            thick_start = start
        if stop < thick_stop:
            thick_stop = stop
        if (start > thick_stop and stop > thick_stop) or (start < thick_start and stop < thick_start):
            thick_start = 0
            thick_stop = 0
        block_starts = ",".join(map(str, block_starts))
        block_sizes = ",".join(map(str, block_sizes))
        return map(str, [self.chromosome, start, stop, name, self.score, self.strand, thick_start, thick_stop,
                         rgb, block_count, block_sizes, block_starts])

    def chromosome_coordinate_to_mrna(self, coord):
        if not (self.start <= coord < self.stop):
            return None
        p = 0
        i = ChromosomeInterval(self.chromosome, coord, coord + 1, self.strand)
        if not any(i.overlap(x) for x in self.exon_intervals):
            return None
        exon_intervals = self.exon_intervals if self.strand == '+' else reversed(self.exon_intervals)
        for e in exon_intervals:
            if i.overlap(e):
                if self.strand == '+':
                    p += coord - e.start
                else:
                    p += e.stop - coord - 1
                break
            p += len(e)
        return p

    def chromosome_coordinate_to_cds(self, coord):
        if not (self.thick_start <= coord < self.thick_stop):
            return None
        p = self.chromosome_coordinate_to_mrna(coord)
        if p is None:
            return p
        return self.mrna_coordinate_to_cds(p)

    def mrna_coordinate_to_chromosome(self, coord):
        if not (0 <= coord < len(self)):
            return None
        p = 0
        exon_intervals = self.exon_intervals if self.strand == '+' else reversed(self.exon_intervals)
        for e in exon_intervals:
            if p + len(e) > coord:
                if self.strand == '+':
                    return e.start + (coord - p)
                else:
                    return e.stop - (coord - p) - 1
            p += len(e)

    def mrna_coordinate_to_cds(self, coord):
        if self.strand == '+':
            cds_start = self.chromosome_coordinate_to_mrna(self.thick_start)
        else:
            cds_start = self.chromosome_coordinate_to_mrna(self.thick_stop - 1)
        r = coord - cds_start
        if not (0 <= r < self.cds_size):
            return None
        return r

    def cds_coordinate_to_mrna(self, coord):
        if not (0 <= coord < self.cds_size):
            return None
        if self.strand == '+':
            cds_start = self.chromosome_coordinate_to_mrna(self.thick_start)
        else:
            cds_start = self.chromosome_coordinate_to_mrna(self.thick_stop - 1)
        return cds_start + coord

    def cds_coordinate_to_chromosome(self, coord):
        if not (0 <= coord < self.cds_size):
            return None
        if self.strand == '+':
            cds_start = self.chromosome_coordinate_to_mrna(self.thick_start)
        else:
            cds_start = self.chromosome_coordinate_to_mrna(self.thick_stop - 1)
        c = self.mrna_coordinate_to_chromosome(cds_start + coord)
        return c

    def get_mrna(self, seq_dict):
        """
        Returns the mRNA sequence for this transcript based on a Fasta object.
        and the start/end positions and the exons. Sequence returned in
        5'-3' transcript orientation.
        """
        sequence = seq_dict[self.chromosome]
        assert self.stop <= len(sequence)
        s = []
        for e in self.exon_intervals:
            s.append(sequence[e.start:e.stop])
        if self.strand == '+':
            mrna = ''.join(s)
        else:
            mrna = reverse_complement(''.join(s))
        return str(mrna)

    def get_sequence(self, seq_dict):
        """
        Returns the entire chromosome sequence for this transcript, (+) strand orientation.
        """
        sequence = seq_dict[self.chromosome]
        return sequence[self.start:self.stop]

    def get_cds(self, seq_dict):
        """
        Return the CDS sequence (as a string) for the transcript
        (based on the exons) using a sequenceDict as the sequence source.
        The returned sequence is in the correct 5'-3' orientation (i.e. it has
        been reverse complemented if necessary).
        """
        sequence = seq_dict[self.chromosome]
        assert self.stop <= len(sequence)
        # make sure this isn't a non-coding gene
        if self.thick_start == self.thick_stop == 0:
            return ''
        s = []
        for e in self.exon_intervals:
            if self.thick_start < e.start and e.stop < self.thick_stop:
                # squarely in the CDS
                s.append(sequence[e.start:e.stop])
            elif e.start <= self.thick_start < e.stop < self.thick_stop:
                # thickStart marks the start of the CDS
                s.append(sequence[self.thick_start:e.stop])
            elif e.start <= self.thick_start and self.thick_stop <= e.stop:
                # thickStart and thickStop mark the whole CDS
                s.append(sequence[self.thick_start: self.thick_stop])
            elif self.thick_start < e.start < self.thick_stop <= e.stop:
                # thickStop marks the end of the CDS
                s.append(sequence[e.start:self.thick_stop])
        if self.strand == '-':
            cds = reverse_complement(''.join(s))
        else:
            cds = ''.join(s)
        return str(cds)

    def get_protein_sequence(self, seq_dict):
        """
        Returns the translated protein sequence for this transcript in single
        character space.
        """
        cds = self.get_cds(seq_dict)
        if len(cds) < 3:
            return ''
        return translate_sequence(self.get_cds(seq_dict).upper())


class GenePredTranscript(Transcript):
    """
    Subclasses Transcript to represent genePred entries. genePred entries have the same information, except that they
    also tell you whether the CDS is complete on both ends, and the frame information of each exon.
    """
    # adding slots for new fields
    __slots__ = ('cds_start_stat', 'cds_end_stat', 'exon_frames', 'name2', 'id')

    def __init__(self, gene_pred_tokens):
        name = gene_pred_tokens[0]
        chrom = gene_pred_tokens[1]
        strand = gene_pred_tokens[2]
        start = gene_pred_tokens[3]
        stop = gene_pred_tokens[4]
        thick_start = gene_pred_tokens[5]
        thick_stop = gene_pred_tokens[6]
        block_count = gene_pred_tokens[7]
        exon_starts = gene_pred_tokens[8]
        exon_ends = gene_pred_tokens[9]
        self.id = gene_pred_tokens[10]
        self.name2 = gene_pred_tokens[11]
        self.cds_start_stat = gene_pred_tokens[12]
        self.cds_end_stat = gene_pred_tokens[13]
        self.exon_frames = [int(x) for x in gene_pred_tokens[14].split(',') if x != '']
        # convert genePred format coordinates to BED-like coordinates to make intervals
        block_starts = [int(x) for x in exon_starts.split(',') if x != '']
        block_ends = [int(x) for x in exon_ends.split(',') if x != '']
        block_sizes = ",".join(map(str, [e - s for e, s in izip(block_ends, block_starts)]))
        block_starts = ",".join(map(str, [x - int(start) for x in block_starts]))
        bed_tokens = [chrom, start, stop, name, '0', strand, thick_start, thick_stop, '0', block_count, block_sizes,
                      block_starts]
        super(GenePredTranscript, self).__init__(bed_tokens)

    @property
    def offset(self):
        frames = [x for x in self.exon_frames if x != -1]
        if len(frames) == 0:
            return 0
        if self.strand == '+':
            offset = 3 - frames[0]
        else:
            offset = 3 - frames[-1]
        if offset == 3:
            offset = 0
        return offset

    def get_protein_sequence(self, seq_dict):
        """
        Returns the translated protein sequence for this transcript in single
        character space.
        """
        cds = self.get_cds(seq_dict, in_frame=True)
        if len(cds) < 3:
            return ""
        try:
            return translate_sequence(cds.upper())
        except AssertionError:
            raise RuntimeError('Failed to translate transcript {} with sequence {}'.format(self.name, cds))

    def get_cds(self, seq_dict, in_frame=False):
        """
        Returns the CDS sequence. Overrides the parental get_cds function to provide frame-corrected sequence.
        Note that if a in-frame sequence is requested, it will no longer correspond with internal coordinates.
        """
        cds = super(GenePredTranscript, self).get_cds(seq_dict)
        if in_frame is False:
            return cds
        else:
            offset = self.offset
            return cds[offset:len(cds) - ((len(cds) - offset) % 3)]

    def get_gene_pred(self, name=None, start_offset=None, stop_offset=None, name2=None, uid=None):
        """
        Returns this transcript as a genePred transcript.
        If start_offset or stop_offset are set (chromosome coordinates), then this record will be changed to only
        show results within that region, which is defined in chromosome coordinates.
        TODO: if the functionality to resize this is used, exon frame information will be incorrect.
        """
        bed_rec = self.get_bed(name=name, new_start=start_offset, new_stop=stop_offset)
        chrom = bed_rec[0]
        start = int(bed_rec[1])
        stop = int(bed_rec[2])
        name = bed_rec[3]
        strand = bed_rec[5]
        thick_start = int(bed_rec[6])
        thick_stop = int(bed_rec[7])
        block_count = int(bed_rec[9])
        block_sizes = map(int, bed_rec[10].split(","))
        block_starts = map(int, bed_rec[11].split(","))
        # convert BED fields to genePred fields
        exon_starts = [start + x for x in block_starts]
        exon_ends = [x + y for x, y in zip(*[exon_starts, block_sizes])]
        exon_starts = ",".join(map(str, exon_starts))
        exon_ends = ",".join(map(str, exon_ends))
        exon_frames = ",".join(map(str, self.exon_frames))
        # change names if desired
        name2 = self.name2 if name2 is None else name2
        uid = self.id if uid is None else uid
        return map(str, [name, chrom, strand, start, stop, thick_start, thick_stop, block_count,
                         exon_starts, exon_ends, uid, name2, self.cds_start_stat, self.cds_end_stat,
                         exon_frames])


def get_gene_pred_dict(gp_file):
    """
    Produces a dictionary of GenePredTranscripts from a genePred file
    :param gp_file: A genePred file path or handle.
    :return: A dictionary of name:transcript pairs
    """
    return {t.name: t for t in gene_pred_iterator(gp_file)}


def gene_pred_iterator(gp_file):
    """
    Iterator for GenePred file or handle, producing tuples of (name, GenePredTranscript)
    :param gp_file: A genePred file path or handle.
    :return: tuples of (name, GenePredTranscript)
    """
    for tokens in iter_lines(gp_file):
        if len(tokens) != 15:
            raise RuntimeError('GenePred line had {} tokens, not 15. Record: {}'.format(len(tokens), tokens))
        t = GenePredTranscript(tokens)
        yield t


def get_transcript_dict(bed_file):
    """
    Produces a dictionary of Transcripts from a BED file
    :param bed_file: A BED file path or handle.
    :return: A dictionary of name:transcript pairs
    """
    return {t.name: t for t in transcript_iterator(bed_file)}


def transcript_iterator(bed_file):
    """
    Iterator for BED file or handle, producing tuples of (name, Transcript)
    :param bed_file: A BED file path or handle.
    :return: tuples of (name, Transcript)
    """
    with open(bed_file) as inf:
        for tokens in iter_lines(inf):
            if len(tokens) != 12:
                raise RuntimeError('BED line had {} tokens, not 12. Record: {}'.format(len(tokens), tokens))
            t = Transcript(tokens)
            yield t


def load_gps(gp_list):
    """helper function that loads a list of genePreds into one mega-dict"""
    r = {}
    for gp in gp_list:
        for t in gene_pred_iterator(gp):
            if t.name in r:
                raise RuntimeError('Attempted to add duplicate GenePredTranscript object with name {}'.format(name))
            r[t.name] = t
    return r
