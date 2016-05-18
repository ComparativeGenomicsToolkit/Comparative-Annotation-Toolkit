"""
Classifiers that evaluate transcripts. Some are single genome, some rely on mapping information.
"""
import tools.transcripts

########################################################################################################################
# Functions and global variables
########################################################################################################################


short_intron_size = 50
short_cds_size = 75


def short_cds(t):
    """
    Many classifiers do not apply to CDS below a cutoff size.
    """
    return True if t.cds_size <= short_cds_size else False


def short_intron(intron):
    """
    Many classifiers rely on analyzing introns either above or below a cutoff size
    """
    return True if len(intron) <= short_intron_size else False


def is_cds(intron, t):
    return intron.start >= t.thick_start and intron.stop < t.thick_stop


def is_not_cds(intron, t):
    return not is_cds(intron, t)


def get_adjusted_starts_ends(t, aln):
    """
    Converts target coordinates to query for each intron start.
    """
    return [[aln.target_coordinate_to_query(intron.start - 1), aln.target_coordinate_to_query(intron.stop)]
            for intron in t.intron_intervals]


def is_fuzzy_intron(intron, aln, ref_starts, fuzz_distance=8):
    """
    Determines if a intron is within fuzz distance of its aligned partner.
    """
    q_gap_start = aln.target_coordinate_to_query(intron.start - 1)
    q_gap_end = aln.target_coordinate_to_query(intron.stop)
    return query_contains_intron(q_gap_start - fuzz_distance, q_gap_end + fuzz_distance, ref_starts)


def query_contains_intron(q_gap_start, q_gap_end, ref_starts):
    r = [q_gap_start <= ref_gap <= q_gap_end for ref_gap in ref_starts]
    return True if any(r) else False


def fix_ref_q_starts(ref_aln):
    """
    Inverts a negative strand reference psl.
    """
    if ref_aln.strand == "-":
        ref_starts = [ref_aln.q_size - (ref_aln.q_starts[i] + ref_aln.block_sizes[i]) for i in
                      xrange(len(ref_aln.q_starts))]
    else:
        ref_starts = ref_aln.q_starts
    return ref_starts


########################################################################################################################
# Single genome
########################################################################################################################


class BadFrame(object):
    """
    Looks for CDS sequences that are not a multiple of 3. Must have at least short_cds_size coding bases.
    Will report a BED record of the transcript if true
    """
    def __call__(self, a, ref_fasta):
        # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
        if short_cds(a):
            return []
        if a.cds_size % 3 != 0:
            return [tools.transcripts.chromosome_coordinate_to_bed(a, a.thick_start, a.thick_stop, self.rgb, self.name)]
        else:
            return []


class BeginStart(object):
    """
    Is the first 3 bases of thick_start 'ATG'?
    Returns a BED record of the first 3 bases if this is NOT true
    """
    def __call__(self, a, ref_fasta):
        # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
        if short_cds(a):
            return []
        elif a.get_cds(ref_fasta)[:3] != "ATG":
            return [tools.transcripts.cds_coordinate_to_bed(a, 0, 3, self.rgb, self.name)]
        else:
            return []


class EndStop(object):
    """
    Are the last three bases a stop codon?
    If this is NOT true, will report a BED record of the last 3 bases.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def __call__(self, a, ref_fasta):
        stop_codons = {'TAA', 'TGA', 'TAG'}
        # do not include noncoding transcripts or lift-overs that contain less than short_cds_size
        if short_cds(a):
            return []
        elif a.get_cds(ref_fasta)[-3:] not in stop_codons:
            return [tools.transcripts.cds_coordinate_to_bed(a, a.cds_size - 3, a.cds_size, self.rgb, self.name)]
        else:
            return []


class CdsGap(object):
    """
    Are any of the CDS introns too short?
    Reports a BED record for each intron interval that is too short.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def __call__(self, a, ref_fasta):
        bed_recs = []
        for intron in a.intron_intervals:
            if len(intron) <= short_intron_size:
                if len(intron) % 3 != 0 and is_cds(intron, a) is True:
                    bed_rec = tools.transcripts.interval_to_bed(a, intron, self.rgb, self.name)
                    bed_recs.append(bed_rec)
        return bed_recs


class CdsMult3Gap(CdsGap):
    """
    Same as CdsGap, but only reports on multiple of 3s.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def __call__(self, a, ref_fasta):
        bed_recs = []
        for intron in a.intron_intervals:
            if short_intron(intron):
                if len(intron) % 3 == 0 and is_cds(intron, a) is True:
                    bed_rec = tools.transcripts.interval_to_bed(a, intron, self.rgb, self.name)
                    bed_recs.append(bed_rec)
        return bed_recs


class UtrGap(CdsGap):
    """
    Are any UTR introns too short?
    Reports on all such introns.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def __call__(self, a, ref_fasta):
        bed_recs = []
        for intron in a.intron_intervals:
            if short_intron(intron) and is_cds(intron, a) is False:
                bed_rec = tools.transcripts.interval_to_bed(a, intron, self.rgb, self.name)
                bed_recs.append(bed_rec)
        return bed_recs


class UnknownGap(CdsGap):
    """
    Looks for short introns that contain unknown bases. Any number of unknown bases is fine.
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def __call__(self, a, ref_fasta):
        bed_recs = []
        for intron in a.intron_intervals:
            if short_intron(intron):
                if "N" in intron.get_sequence(ref_fasta):
                    bed_rec = tools.transcripts.interval_to_bed(a, intron, self.rgb, self.name)
                    bed_recs.append(bed_rec)
        return bed_recs


class CdsNonCanonSplice(object):
    """
    Are any of the CDS introns splice sites not of the canonical form
    GT..AG
    Ignores cases where the splice sites are ambiguous (contains an N)
    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def __call__(self, a, ref_fasta, cds_filter_fn=is_cds, splice_dict={"GT": "AG"}):
        bed_recs = []
        for intron in a.intron_intervals:
            splice_is_good = analyze_splice(intron, a, ref_fasta, cds_filter_fn, splice_dict)
            if splice_is_good is True:
                bed_rec = tools.transcripts.splice_intron_interval_to_bed(a, intron, self.rgb, self.name)
                bed_recs.append(bed_rec)
        return bed_recs


class CdsUnknownSplice(CdsNonCanonSplice):
    """
    Are any of the CDS introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    def __call__(self, a, ref_fasta, cds_filter_fn=is_cds, splice_dict={"GT": "AG", "GC": "AG", "AT": "AC"}):
        return CdsNonCanonSplice.__call__(self, a, ref_fasta, cds_filter_fn, splice_dict)


class UtrNonCanonSplice(CdsNonCanonSplice):
    """
    Are any of the UTR introns splice sites not of the canonical form
    GT..AG

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    def __call__(self, a, ref_fasta, cds_filter_fn=is_not_cds, splice_dict={"GT": "AG"}):
        return CdsNonCanonSplice.__call__(self, a, ref_fasta, cds_filter_fn, splice_dict)


class UtrUnknownSplice(CdsNonCanonSplice):
    """
    Are any of the UTR introns splice sites not of the form
    GT..AG, GC..AG, AT..AC

    This classifier is only applied to introns which are longer than
    a minimum intron size.
    """
    def __call__(self, a, ref_fasta, cds_filter_fn=is_not_cds, splice_dict={"GT": "AG", "GC": "AG", "AT": "AC"}):
        return CdsNonCanonSplice.__call__(self, a, ref_fasta, cds_filter_fn, splice_dict)


class SpliceContainsUnknownBases(object):
    """
    Do any of the splice junctions contain unknown bases?
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def __call__(self, a, ref_fasta):
        bed_recs = []
        for intron in a.intron_intervals:
            if short_intron(intron) is False:
                seq = intron.get_sequence(ref_fasta, stranded=True)
                donor, acceptor = seq[:2], seq[-2:]
                if "N" in donor or "N" in acceptor:
                    bed_rec = tools.transcripts.splice_intron_interval_to_bed(a, intron, self.rgb, self.name)
                    bed_recs.append(bed_rec)
        return bed_recs


class InFrameStop(object):
    """
    Reports on in frame stop codons for each transcript.

    In order to be considered, must have at least 25 codons.

    Returns a BED record of the position of the in frame stop if it exists.
    """
    @property
    def rgb(self):
        return self.colors["mutation"]

    def __call__(self, a, ref_fasta):
        bed_recs = []
        cds = a.get_cds(ref_fasta, in_frame=False)
        offset = tools.transcripts.find_offset(a.exon_frames, a.strand)
        for i, codon in bio_lib.read_codons_with_position(cds, offset, skip_last=True):
            amino_acid = bio_lib.codon_to_amino_acid(codon)
            if amino_acid == "*":
                bed_rec = tools.transcripts.cds_coordinate_to_bed(a, i, i + 3, self.rgb, self.name)
                bed_recs.append(bed_rec)
        return bed_recs


class ShortCds(object):
    """
    Looks to see if this transcript has a short CDS.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def __call__(self, a, ref_fasta):
        if short_cds(a) is True and a.cds_size != 0:
            bed_rec = tools.transcripts.cds_coordinate_to_bed(a, 0, a.cds_size, self.rgb, self.name)
            return [bed_rec]
        else:
            return []


class UnknownBases(object):
    """
    Does this alignment contain Ns in the target genome?

    Only looks at mRNA bases, and restricts to CDS if cds is True
    """
    @property
    def rgb(self):
        return self.colors["assembly"]

    def make_bed_recs(self, a, s, bed_rec_fn, r=re.compile("[atgcATGC][N]+[atgcATGC]")):
        for m in re.finditer(r, s):
            yield bed_rec_fn(a, m.start() + 1, m.end() - 1, self.rgb, self.name)

    def __call__(self, a, ref_fasta, cds=False):
        if cds is True:
            bed_rec_fn = tools.transcripts.cds_coordinate_to_bed
            s = a.get_cds(ref_fasta)
        else:
            bed_rec_fn = tools.transcripts.transcript_coordinate_to_bed
            s = a.get_mrna(ref_fasta)
        bed_recs = list(self.make_bed_recs(a, s, bed_rec_fn))
        return bed_recs


class UnknownCdsBases(UnknownBases):
    def __call__(self, a, ref_fasta, cds=True):
        return UnknownBases.__call__(self, a, ref_fasta, cds)


class LongTranscript(object):
    """
    Is this transcript unbelievably long? Filters out poor alignments.
    """
    @property
    def rgb(self):
        return self.colors["alignment"]

    def __call__(self, a, ref_fasta, size_cutoff=3 * 10 ** 6):
        if a.stop - a.start >= size_cutoff:
            return [tools.transcripts.transcript_to_bed(a, self.rgb, self.name)]
        else:
            return []
