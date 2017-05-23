"""
Represent either BED12 or genePred transcripts as objects. Allows for conversion of coordinates between
chromosome, mRNA and CDS coordinate spaces. Can slice objects into subsets.
"""
import collections
from itertools import izip
from bx.intervals.cluster import ClusterTree

from mathOps import find_closest, find_intervals
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

    def __repr__(self):
        return 'Transcript({})'.format(self.get_bed())

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

    def get_gtf_coding_interval(self):
        """
        GTF/GFF3 files do not include the stop codon as part of the CDS. For this reason, we have a special case.

        Using the get_bed() functionality to achieve this.
        """
        if self.strand == '+':
            # special case -- we have a 3bp transcript with a complete start marking. Then, start codon == stop codon
            if self.cds_size == 3:
                bed = self.get_bed(new_stop=self.cds_coordinate_to_chromosome(3), new_start=self.thick_start)
            else:
                bed = self.get_bed(new_stop=self.cds_coordinate_to_chromosome(self.cds_size - 3),
                                   new_start=self.thick_start)
        else:
            bed = self.get_bed(new_start=self.cds_coordinate_to_chromosome(self.cds_size - 3),
                               new_stop=self.thick_stop)
        return Transcript(bed).coding_interval

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
        else:
            new_start = self.start
        if new_stop is not None:
            assert new_stop <= self.stop
        else:
            new_stop = self.stop
        rgb = self.rgb if rgb is None else rgb
        name = self.name if name is None else name

        # special case -- start == stop
        if new_start == new_stop:
            if self.cds_size == 0:
                thick_start = thick_stop = 0
            else:
                thick_start = new_start
                thick_stop = new_stop
            return map(str, [self.chromosome, new_start, new_stop, name, self.score, self.strand, thick_start,
                             thick_stop, rgb, 1, 0, 0])

        if self.chromosome_coordinate_to_mrna(new_start) is None:
            new_start = find_closest([x.start for x in self.exon_intervals], new_start)
        if self.chromosome_coordinate_to_mrna(new_stop) is None:
            new_stop = find_closest([x.stop for x in self.exon_intervals], new_stop)

        # start slicing out intervals
        new_interval = ChromosomeInterval(self.chromosome, new_start, new_stop, self.strand)
        exon_intervals = []
        for exon in self.exon_intervals:
            new_exon = exon.intersection(new_interval)
            if new_exon is None:
                continue
            exon_intervals.append(new_exon)

        # if new_start or new_stop were not within the exonic intervals, adjust them
        if new_start != exon_intervals[0].start:
            new_start = exon_intervals[0].start
        if new_stop != exon_intervals[-1].stop:
            new_stop = exon_intervals[-1].stop
        thick_start = max(self.thick_start, new_start)
        thick_stop = min(self.thick_stop, new_stop)
        if thick_start >= self.thick_stop or thick_stop < self.thick_start:
            thick_start = 0
            thick_stop = 0
        block_count = len(exon_intervals)
        block_sizes = ','.join(map(str, [len(x) for x in exon_intervals]))
        block_starts = ','.join(map(str, [x.start - new_start for x in exon_intervals]))
        return map(str, [self.chromosome, new_start, new_stop, name, self.score, self.strand, thick_start, thick_stop,
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

    def get_start_intervals(self):
        """
        Returns one or more ChromosomeInterval objects that represents the starting CDS interval for this transcript.
        More than one may exist if the codon is split over a splice junction.
        """
        assert self.cds_size >= 3
        positions = sorted([self.cds_coordinate_to_chromosome(x) for x in range(3)])
        merged_intervals = list(find_intervals(positions))
        intervals = [ChromosomeInterval(self.chromosome, i[0], i[-1] + 1, self.strand) for i in merged_intervals]
        assert sum(len(x) for x in intervals) == 3
        c = 0
        for i in intervals:
            i.data = convert_frame(c)
            c += len(i)
        return intervals

    def get_stop_intervals(self):
        """
        Returns one or more ChromosomeInterval objects that represents the ending CDS interval for this transcript.
        More than one may exist if the codon is split over a splice junction.
        """
        assert self.cds_size >= 3
        positions = sorted([self.cds_coordinate_to_chromosome(x) for x in range(self.cds_size - 3, self.cds_size)])
        merged_intervals = list(find_intervals(positions))
        intervals = [ChromosomeInterval(self.chromosome, i[0], i[-1] + 1, self.strand) for i in merged_intervals]
        assert sum(len(x) for x in intervals) == 3
        c = 0
        for i in intervals:
            i.data = convert_frame(c)
            c += len(i)
        return intervals


class GenePredTranscript(Transcript):
    """
    Subclasses Transcript to represent genePred entries. genePred entries have the same information, except that they
    also tell you whether the CDS is complete on both ends, and the frame information of each exon.
    """
    # adding slots for new fields
    __slots__ = ('cds_start_stat', 'cds_end_stat', 'exon_frames', 'name2', 'score')

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
        self.score = gene_pred_tokens[10]
        self.name2 = gene_pred_tokens[11]
        self.cds_start_stat = gene_pred_tokens[12]
        self.cds_end_stat = gene_pred_tokens[13]
        self.exon_frames = [int(x) for x in gene_pred_tokens[14].split(',') if x != '']
        # convert genePred format coordinates to BED-like coordinates to make intervals
        block_starts = [int(x) for x in exon_starts.split(',') if x != '']
        block_ends = [int(x) for x in exon_ends.split(',') if x != '']
        block_sizes = ",".join(map(str, [e - s for e, s in izip(block_ends, block_starts)]))
        block_starts = ",".join(map(str, [x - int(start) for x in block_starts]))
        bed_tokens = [chrom, start, stop, name, self.score, strand, thick_start, thick_stop, '0', block_count,
                      block_sizes, block_starts]
        super(GenePredTranscript, self).__init__(bed_tokens)

    def __repr__(self):
        return 'GenePredTranscript({})'.format(self.get_gene_pred())

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

    def _get_exon_intervals(self):
        """
        Overrides _get_exon_intervals to attach frame information to the intervals
        :return: List of ChromosomeIntervals
        """
        exon_intervals = []
        for block_size, block_start, frame in izip(*(self.block_sizes, self.block_starts, self.exon_frames)):
            start = self.start + block_start
            stop = self.start + block_start + block_size
            exon_intervals.append(ChromosomeInterval(self.chromosome, start, stop, self.strand, data={'frame': frame}))
        return exon_intervals

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

    def get_gene_pred(self, name=None, new_start=None, new_stop=None, name2=None, score=None):
        """
        Returns this transcript as a genePred transcript.
        If new_start or new_stop are set (chromosome coordinates), then this record will be changed to only
        show results within that region, which is defined in chromosome coordinates. The frames field will be properly
        adjusted, and the cds_start_stat/cds_end_stat fields will change to 'unk' if they are moved

        TODO: If this is a transMap transcript, and there were coding indels, the frame information will change to
        reflect the new arrangement and the implicit indel information will be lost.
        """
        name = self.name if name is None else name
        name2 = self.name2 if name2 is None else name2
        score = self.score if score is None else score

        # if no resizing, just return what we have
        if new_start is None and new_stop is None:
            exon_starts = ','.join(map(str, [exon.start for exon in self.exon_intervals]))
            exon_ends = ','.join(map(str, [exon.stop for exon in self.exon_intervals]))
            exon_frames = ','.join(map(str, self.exon_frames))
            return map(str, [name, self.chromosome, self.strand, self.start, self.stop, self.thick_start,
                             self.thick_stop, len(self.exon_intervals), exon_starts, exon_ends, score, name2,
                             self.cds_start_stat, self.cds_end_stat, exon_frames])
        if new_start is not None and new_stop is not None:
            assert new_start <= new_stop
        if new_start is not None:
            assert new_start >= self.start
        else:
            new_start = self.start
        if new_stop is not None:
            assert new_stop <= self.stop
        else:
            new_stop = self.stop

        # start slicing out intervals, adjusting the frames
        new_interval = ChromosomeInterval(self.chromosome, new_start, new_stop, self.strand)
        exon_intervals = []
        exon_frames = []
        exon_iter = self.exon_intervals if self.strand == '+' else self.exon_intervals[::-1]
        frame_iter = self.exon_frames if self.strand == '+' else reversed(self.exon_frames)

        # attempt to find the first frame. If there is none, then we have a non-coding transcript and this is easy
        try:
            starting_frame = [f for f in frame_iter if f != -1][0]
        except IndexError:  # non-coding transcript
            exon_intervals = [exon.intersection(new_interval) for exon in exon_iter]
            exon_frames = [-1] * len(exon_intervals)
        else:  # start following frame to adjust for resized transcript
            cds_counter = 0  # keep track of total CDS bases encountered
            cds_flag = False
            for exon in exon_iter:
                new_exon = exon.intersection(new_interval)
                if new_exon is None:
                    continue
                exon_intervals.append(new_exon)
                coding_exon = exon.intersection(self.coding_interval)
                if coding_exon is None:
                    exon_frames.append(-1)
                elif cds_flag is False:
                    cds_flag = True
                    exon_frames.append(starting_frame)
                    cds_counter += len(coding_exon) + starting_frame
                else:
                    exon_frames.append(cds_counter % 3)
                    cds_counter += len(coding_exon)

        # flip back around negative strand transcripts
        if self.strand == '-':
            exon_intervals = exon_intervals[::-1]
            exon_frames = exon_frames[::-1]

        # if new_start or new_stop were intronic coordinates, fix this
        if new_start != exon_intervals[0].start:
            new_start = exon_intervals[0].start
        if new_stop != exon_intervals[-1].stop:
            new_stop = exon_intervals[-1].stop

        thick_start = max(self.thick_start, new_start)
        thick_stop = min(self.thick_stop, new_stop)
        cds_start_stat = 'unk' if thick_start != self.thick_start else self.cds_start_stat
        cds_end_stat = 'unk' if thick_stop != self.thick_stop else self.cds_end_stat
        exon_count = len(exon_intervals)
        exon_starts = ','.join(map(str, [exon.start for exon in exon_intervals]))
        exon_ends = ','.join(map(str, [exon.stop for exon in exon_intervals]))
        exon_frames = ','.join(map(str, exon_frames))
        return map(str, [name, self.chromosome, self.strand, new_start, new_stop, thick_start, thick_stop, exon_count,
                         exon_starts, exon_ends, score, name2, cds_start_stat, cds_end_stat, exon_frames])


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
                raise RuntimeError('Attempted to add duplicate GenePredTranscript object with name {}'.format(t.name))
            r[t.name] = t
    return r


def convert_frame(exon_frame):
    """converts genePred-style exonFrame to GFF-style phase"""
    mapping = {0: 0, 1: 2, 2: 1, -1: '.'}
    return mapping[exon_frame]


def create_bed_info_gp(gp):
    """Creates the block_starts, block_sizes and exon_frames fields from a GenePredTranscript object"""
    block_starts = ','.join(map(str, gp.block_starts))
    block_sizes = ','.join(map(str, gp.block_sizes))
    exon_frames = ','.join(map(str, gp.exon_frames))
    return block_starts, block_sizes, exon_frames


def group_transcripts_by_name2(tx_iter):
    """Takes a iterable of GenePredTranscript objects and groups them by name2"""
    r = collections.defaultdict(list)
    for tx in tx_iter:
        r[tx.name2].append(tx)
    return r


def intervals_to_bed(intervals, name=None, score=0, rgb=0, thick_start=0, thick_stop=0):
    """Converts an iterable of intervals into a Transcript object. If any intervals overlap this will fail"""
    assert len(set(i.strand for i in intervals)) == 1
    assert len(set(i.chromosome for i in intervals)) == 1
    intervals = sorted(intervals)
    start = intervals[0].start
    stop = intervals[-1].stop
    block_sizes = ','.join(map(str, [len(i) for i in intervals]))
    block_starts = ','.join(map(str, [i.start - start for i in intervals]))
    i = intervals[0]
    return Transcript([i.chromosome, start, stop, name, score, i.strand, thick_start, thick_stop, rgb,
                       len(intervals), block_sizes, block_starts])


def cluster_txs(txs):
    """Uses a ClusterTree to cluster to cluster transcript objects"""
    cluster_trees = collections.defaultdict(lambda: ClusterTree(0, 1))
    for i, tx in enumerate(txs):
        cluster_trees[tx.chromosome].insert(tx.start, tx.stop, i)
    # convert the clusters to a nested structure of chrom -> cluster_id -> tx objects
    clustered_reads = collections.defaultdict(dict)
    cluster_id = 0
    for chrom, cluster_tree in cluster_trees.iteritems():
        for start, end, interval_indices in cluster_tree.getregions():
            clustered_reads[chrom][cluster_id] = [txs[ix] for ix in interval_indices]
            cluster_id += 1
    return clustered_reads


def divide_clusters(clustered_reads, ref_names):
    """
    Takes the output of cluster_txs and splits them into two groups based on having their name be in ref_names or not.

    Returns a dict mapping cluster IDs to tuples of [ref_txs, non_ref_txs].

    Discards any cluster that does not contain members of both ref and non-ref.

    """
    divided_clusters = {}
    for chrom in clustered_reads:
        for cluster_id, tx_list in clustered_reads[chrom].iteritems():
            ref = [tx for tx in tx_list if tx.name in ref_names and len(tx.intron_intervals) > 0]
            iso = [tx for tx in tx_list if tx.name not in ref_names and len(tx.intron_intervals) > 0]
            if len(ref) > 0 and len(iso) > 0:
                divided_clusters[cluster_id] = [ref, iso]
    return divided_clusters


def construct_start_stop_intervals(intron_intervals, d):
    """Splits a iterable of intervals into two parallel tuples of 0bp intervals representing their start and stop"""
    left_intervals = []
    right_intervals = []
    for i in intron_intervals:
        left_intervals.append(ChromosomeInterval(i.chromosome, i.start - d, i.start + d, i.strand))
        right_intervals.append(ChromosomeInterval(i.chromosome, i.stop - d, i.stop + d, i.strand))
    return tuple(left_intervals), tuple(right_intervals)


def find_subset_match(iso_intervals, enst_intervals):
    """
    Compares intervals produced by construct_start_stop_intervals to each other to find subset matches.
    Used for fuzzy matching of IsoSeq transcripts (iso_intervals) to existing annotations (enst_intervals)
    """
    iso_l, iso_r = iso_intervals
    enst_l, enst_r = enst_intervals
    # if we have fewer reference junctions than isoseq, we can't have a subset match by definition
    if len(iso_l) > len(enst_l):
        return False
    lm = all([any([il.overlap(el) for el in enst_l]) for il in iso_l])
    lr = all([any([ir.overlap(er) for er in enst_r]) for ir in iso_r])
    return lm and lr


def calculate_subset_matches(divided_clusters, fuzz_distance=8):
    """
    A wrapper for find_subset_match that looks at every cluster of transcripts produced by divide_clusters and finds
    a fuzzy match between any non-reference sequence and a reference sequence.

    """
    r = collections.defaultdict(list)
    for cluster_id, (ensts, isos) in divided_clusters.iteritems():
        enst_intervals = collections.defaultdict(list)
        for tx in ensts:
            enst_interval = construct_start_stop_intervals(tx.intron_intervals, fuzz_distance)
            enst_intervals[tuple(enst_interval)].append(tx)
        for iso in isos:
            iso_intervals = construct_start_stop_intervals(iso.intron_intervals, fuzz_distance)
            for enst_interval, enst_txs in enst_intervals.iteritems():
                m = find_subset_match(iso_intervals, enst_interval)
                if m:
                    r[iso.name].extend(enst_txs)
    return r
