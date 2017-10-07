"""
A series of classifiers that evaluate transMap, AugustusTMR and AugustusCGP output.

These classifiers are broken down into 2 groups, which will each end up as a table in the database:

<alnMode>_<txMode>_Metrics:

These classifiers are per-transcript evaluations based on both the transcript alignment and the genome context.
1. PercentUnknownBases: % of mRNA bases that are Ns.
2. AlnCoverage: Alignment coverage in transcript space.
3. AlnIdentity: Alignment identity in transcript space.
4. OriginalIntronVector: Number of original introns within a wiggle distance of any introns in the target. Different
   from the intron support calculated by homGeneMapping because it verifies the actual isoform.
5. ValidStart -- start with ATG?
6. ValidStop -- valid stop codon (in frame)?
7. ProperOrf -- is the orf a multiple of 3?

<alnMode>_<txMode>_Evaluation:

These classifiers are per-transcript evaluations based on the transcript alignment.
Unlike the other two tables, this table stores the actual location of the problems (in genome coordinates) as a
BED-like format. In cases where there are multiple problems, they will be additional rows.
1. CodingInsertion: Do we have any frame-shifting coding insertions?
2. CodingDeletion: Do we have any frame-shifting coding deletions?
3. CodingMult3Insertion: Do we have any mod3 coding insertions?
4. CodingMult3Deletion: Do we have any mod3 coding deletions?
5. NonCodingInsertion: Do we have indels in UTR sequence?
6. NonCodingDeletion: Do we have any indels in UTR sequence?
7. InFrameStop: Are there any in-frame stop codons?


The Metrics and Evaluation groups will have multiple tables for each of the input methods used:
txMode:
1) transMap
2) augTM
3) augTMR
4) augCGP
5) augPB

alnMode:
1) CDS
2) mRNA

"""
import itertools

import pandas as pd

import tools.bio
import tools.dataOps
import tools.fileOps
import tools.intervals
import tools.mathOps
import tools.nameConversions
import tools.psl
import tools.sqlInterface
import tools.transcripts

# hard coded variables
# fuzz distance is the distance between introns allowed in intron coordinates before triggering NumMissingIntrons
# fuzz distance is counted from both sides of the intron
fuzz_distance = 7


def classify(eval_args):
    """
    Runs alignment classification for all alignment modes
    :param eval_args: argparse Namespace produced by EvaluateTranscripts.get_args()
    :return: list of (tablename, dataframe) tuples
    """
    # load shared inputs
    ref_tx_dict = tools.transcripts.get_gene_pred_dict(eval_args.annotation_gp)
    tx_biotype_map = tools.sqlInterface.get_transcript_biotype_map(eval_args.ref_db_path)
    seq_dict = tools.bio.get_sequence_dict(eval_args.fasta)
    # results stores the final dataframes
    results = {}
    for tx_mode, path_dict in eval_args.transcript_modes.iteritems():
        tx_dict = tools.transcripts.get_gene_pred_dict(path_dict['gp'])
        aln_modes = ['CDS', 'mRNA'] if tx_mode != 'augCGP' else ['CDS']
        for aln_mode in aln_modes:
            psl_iter = list(tools.psl.psl_iterator(path_dict[aln_mode]))
            mc_df = metrics_classify(aln_mode, ref_tx_dict, tx_dict, tx_biotype_map, psl_iter, seq_dict)
            ec_df = evaluation_classify(aln_mode, ref_tx_dict, tx_dict, tx_biotype_map, psl_iter, seq_dict)
            results[tools.sqlInterface.tables[aln_mode][tx_mode]['metrics'].__tablename__] = mc_df
            results[tools.sqlInterface.tables[aln_mode][tx_mode]['evaluation'].__tablename__] = ec_df
    return results


def metrics_classify(aln_mode, ref_tx_dict, tx_dict, tx_biotype_map, psl_iter, seq_dict):
    """
    Calculates the alignment metrics and the number of missing original introns on this transcript_chunk
    :return: DataFrame
    """
    r = []
    for ref_tx, tx, psl, biotype in tx_iter(psl_iter, ref_tx_dict, tx_dict, tx_biotype_map):
        original_intron_vector = calculate_original_intron_vector(ref_tx, tx, psl, aln_mode)
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'AlnCoverage', 100 * psl.target_coverage])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'AlnIdentity', 100 * psl.identity])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'AlnGoodness', 100 * (1 - psl.badness)])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'PercentUnknownBases', psl.percent_n])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'OriginalIntrons', original_intron_vector])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'ValidStart', tools.transcripts.has_start_codon(seq_dict, tx)])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'ValidStop', tools.transcripts.has_stop_codon(seq_dict, tx)])
        r.append([ref_tx.name2, ref_tx.name, tx.name, 'ProperOrf', tx.cds_size % 3 == 0])
    columns = ['GeneId', 'TranscriptId', 'AlignmentId', 'classifier', 'value']
    df = pd.DataFrame(r, columns=columns)
    df = df.sort_values(columns)
    df = df.set_index('AlignmentId')
    assert len(r) == len(df)
    return df


def evaluation_classify(aln_mode, ref_tx_dict, tx_dict, tx_biotype_map, psl_iter, seq_dict):
    """
    Calculates the evaluation metrics on this transcript_chunk
    :return: DataFrame
    """
    r = []
    for ref_tx, tx, psl, biotype in tx_iter(psl_iter, ref_tx_dict, tx_dict, tx_biotype_map):
        r.extend(find_indels(tx, psl, aln_mode))
        if biotype == 'protein_coding':
            line = in_frame_stop(tx, seq_dict)
            if line is not None:
                r.append(line)
    columns = ['AlignmentId', 'chromosome', 'start', 'stop', 'name', 'score', 'strand', 'thickStart',
               'thickStop', 'rgb', 'blockCount', 'blockSizes', 'blockStarts']
    df = pd.DataFrame(r, columns=columns)
    df = df.sort_values(columns)
    df = df.set_index('AlignmentId')
    assert len(r) == len(df)
    return df


###
# Metrics Classifiers
###


def calculate_original_intron_vector(ref_tx, tx, psl, aln_mode):
    """
    Determines how many of the gaps present in a given transcript are within a wiggle distance of the parent.

    Algorithm:
    1) Convert the coordinates of each block in the transcript to mRNA/CDS depending on the alignment.
    2) Use the mRNA/CDS alignment to calculate a mapping between alignment positions and transcript positions.
    3) Determine if each block gap coordinate is within fuzz_distance of a parental block gap.

    :param ref_tx: GenePredTranscript object representing the parent transcript
    :param tx: GenePredTranscript object representing the target transcript
    :param psl: PslRow object representing the mRNA/CDS alignment between ref_tx and tx
    :param aln_mode: One of ('CDS', 'mRNA'). Determines if we aligned CDS or mRNA.
    :return: list
    """
    # before we calculate anything, make sure we have introns to lose
    if len(ref_tx.intron_intervals) == 0:
        return None

    # generate a list of reference introns in current coordinates (mRNA or CDS)
    ref_introns = get_intron_coordinates(ref_tx, aln_mode)

    # generate a list of target introns in current coordinates (mRNA or CDS)
    # note that since this PSL is target-referenced, we use query_coordinate_to_target()
    tgt_introns = []
    for intron in get_intron_coordinates(tx, aln_mode):
        p = psl.query_coordinate_to_target(intron)
        if p is not None:
            tgt_introns.append(p)

    # if we lost all introns due to CDS filtering, return a vector of all 0s
    if len(tgt_introns) == 0:
        return ','.join(['0'] * len(ref_tx.intron_intervals))

    # count the number of introns within wiggle distance of each other
    intron_vector = []
    for ref_intron in ref_introns:
        closest = tools.mathOps.find_closest(tgt_introns, ref_intron)
        if closest - fuzz_distance < ref_intron < closest + fuzz_distance:
            intron_vector.append(1)
        else:
            intron_vector.append(0)
    return ','.join(map(str, intron_vector))


###
# Alignment Evaluation Classifiers
###


def in_frame_stop(tx, fasta):
    """
    Finds the first in frame stop of this transcript, if there are any

    :param tx: Target GenePredTranscript object
    :param fasta: pyfasta Fasta object mapping the genome fasta for this analysis
    :return: A ChromosomeInterval object if an in frame stop was found otherwise None
    """
    for start_pos, stop_pos, codon in tx.codon_iterator(fasta):
        if tools.bio.translate_sequence(codon) == '*':
            bed = tx.get_bed(new_start=start_pos, new_stop=stop_pos, rgb='135,78,191', name='InFrameStop')
            return [tx.name] + bed


def find_indels(tx, psl, aln_mode):
    """
    Walks the psl alignment looking for alignment gaps. Reports all such gaps in Chromosome Coordinates, marking
    the type of gap (CodingInsertion, CodingMult3Insertion, CodingDeletion, CodingMult3Deletion)

    Insertion/Deletion is relative to the target genome, for example:

    CodingInsertion:
    ref: ATGC--ATGC
    tgt: ATGCGGATGC

    CodingDeletion:
    ref: ATGCGGATGC
    tgt: ATGC--ATGC

    :param tx: GenePredTranscript object representing the target transcript
    :param psl: PslRow object describing CDS alignment between ref_tx and tx.
    :param aln_mode: One of ('CDS', 'mRNA'). Determines if we aligned CDS or mRNA.
    :return: list of bed12-format lists
    """
    def convert_coordinates_to_chromosome(left_pos, right_pos, coordinate_fn, strand):
        """convert alignment coordinates to target chromosome coordinates, inverting if negative strand"""
        left_chrom_pos = coordinate_fn(left_pos)
        assert left_chrom_pos is not None
        right_chrom_pos = coordinate_fn(right_pos)
        assert right_chrom_pos is not None
        if strand == '-':
            left_chrom_pos, right_chrom_pos = right_chrom_pos, left_chrom_pos
        assert right_chrom_pos >= left_chrom_pos
        return left_chrom_pos, right_chrom_pos

    def parse_indel(left_pos, right_pos, coordinate_fn, tx, offset, gap_type):
        """Converts either an insertion or a deletion into a output transcript"""
        left_chrom_pos, right_chrom_pos = convert_coordinates_to_chromosome(left_pos, right_pos, coordinate_fn,
                                                                            tx.strand)
        if left_chrom_pos is None or right_chrom_pos is None:
            assert aln_mode == 'CDS'
            return None

        if left_chrom_pos > tx.thick_start and right_chrom_pos < tx.thick_stop:
            indel_type = 'CodingMult3' if offset % 3 == 0 else 'Coding'
        else:
            indel_type = 'NonCoding'

        new_bed = tx.get_bed(new_start=left_chrom_pos, new_stop=right_chrom_pos, rgb=offset,
                             name=''.join([indel_type, gap_type]))
        return [tx.name] + new_bed

    # depending on mode, we convert the coordinates from either CDS or mRNA
    # we also have a different position cutoff to make sure we are not evaluating terminal gaps
    if aln_mode == 'CDS':
        coordinate_fn = tx.cds_coordinate_to_chromosome
    else:
        coordinate_fn = tx.mrna_coordinate_to_chromosome

    # r holds the output
    r = []

    # remember where we were last iteration
    q_pos = 0
    t_pos = 0
    # iterate over block starts[i], q_starts[i + 1], t_starts[i + 1]
    for block_size, q_start, t_start in itertools.izip(*[psl.block_sizes, psl.q_starts[1:], psl.t_starts[1:]]):
        q_offset = q_start - block_size - q_pos
        t_offset = t_start - block_size - t_pos
        assert not (q_offset == t_offset == 0)
        assert (q_offset >= 0 and t_offset >= 0)
        if q_offset != 0:  # query insertion -> insertion in target sequence
            left_pos = q_start - q_offset
            right_pos = q_start
            row = parse_indel(left_pos, right_pos, coordinate_fn, tx, q_offset, 'Insertion')
            if row is not None:
                r.append(row)
        if t_offset != 0:  # target insertion -> insertion in reference sequence
            left_pos = right_pos = q_start
            row = parse_indel(left_pos, right_pos, coordinate_fn, tx, t_offset, 'Deletion')
            if row is not None:
                r.append(row)
        q_pos = q_start
        t_pos = t_start
    return r


###
# Helper functions
###


def tx_iter(psl_iter, ref_tx_dict, tx_dict, tx_biotype_map):
    """
    yields tuples of (GenePredTranscript <reference> , GenePredTranscript <target>, PslRow, biotype
    """
    for psl in psl_iter:
        # this psl is target-referenced
        ref_tx = ref_tx_dict[psl.t_name]
        tx = tx_dict[psl.q_name]
        biotype = tx_biotype_map[psl.t_name]
        yield ref_tx, tx, psl, biotype


def convert_cds_frames(ref_tx, tx, aln_mode):
    """
    Wrapper for convert_cds_frame that converts the reference and target GenePredTranscript objects to CDS-frame
    Transcript objects only if the biotype is protein_coding and the transcripts are out of frame

    :param ref_tx: Reference GenePredTranscript object
    :param tx: Target GenePredTranscript object
    :param aln_mode: If we are in CDS mode, we need to convert the transcripts to a CDS-framed object.
    :return: tuple of GenePredTranscript objects (ref_tx, tx)
    """
    if aln_mode == 'CDS':
        if ref_tx.offset != 0:
            ref_tx = convert_cds_frame(ref_tx)
        if tx.offset != 0:
            tx = convert_cds_frame(tx)
    return ref_tx, tx


def convert_cds_frame(tx):
    """
    If this GenePredTranscript object is out of frame, return a new Transcript object representing just the CDS, in
    frame, trimmed to be a multiple of 3

    :param tx: GenePredTranscript object
    :return: Transcript object
    """
    offset = tx.offset
    mod3 = (tx.cds_size - offset) % 3
    if tx.strand == '+':
        b = tx.get_bed(new_start=tx.thick_start + offset, new_stop=tx.thick_stop - mod3)
    else:
        b = tx.get_bed(new_start=tx.thick_start + mod3, new_stop=tx.thick_stop - offset)
    return tools.transcripts.Transcript(b)


def get_intron_coordinates(tx, aln_mode):
    """
    Converts the block_starts coordinates to mRNA or CDS coordinates used in the alignment based on the alignment mode.

    :param tx:GenePredTranscript object
    :param aln_mode: One of ('CDS', 'mRNA'). Used to determine if we aligned in CDS space or mRNA space
    :return: list of integers
    """
    if aln_mode == 'CDS':
        tx = convert_cds_frame(tx)
        introns = [tx.chromosome_coordinate_to_cds(tx.start + x) for x in tx.block_starts[1:]]
    else:
        introns = [tx.chromosome_coordinate_to_mrna(tx.start + x) for x in tx.block_starts[1:]]
    # remove None which means this transcript is protein_coding and that exon is entirely non-coding
    return [x for x in introns if x is not None]


def get_exon_intervals(tx, aln_mode):
    """
    Generates a dict of intervals for this transcript in either mRNA coordinates or CDS coordinates depending on
    alignment mode.

    We maintain a mapping of where the exon came from to deal with CDS conversions and negative strands easily.

    :param tx: GenePredTranscript object
    :param aln_mode: One of ('CDS', 'mRNA'). Used to determine if we aligned in CDS space or mRNA space
    :return: dict of ChromosomeInterval objects {reference:converted}
    """
    if aln_mode == 'CDS':
        tx = convert_cds_frame(tx)
    exons = {}
    for exon in tx.exon_intervals:
        start = tx.chromosome_coordinate_to_mrna(exon.start)
        stop = tx.chromosome_coordinate_to_mrna(exon.stop - 1)  # zero based, half open
        if tx.strand == '-':
            start, stop = stop, start
        i = tools.intervals.ChromosomeInterval(None, start, stop + 1, '.')
        exons[exon] = i
    return exons
