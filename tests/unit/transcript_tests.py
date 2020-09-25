import unittest
from tools.transcripts import Transcript, GenePredTranscript


class PositiveStrandTranscriptTests(unittest.TestCase):
    """
    Tests the Transcript functionality part of sequence_lib.
    Tests the example positive strand BED record drawn out below:
    chrom    0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
    seq      G  T  A  T  T  C  T  T  G  G  A  C  C  T  A  A
    tx       -  -  a  t  T  C  -  T  G  G  -  -  C  t  a  -
    tx.pos         0  1  2  3     4  5  6        7  8  9
    cds.pos              0  1     2  3  4        5
    """

    def setUp(self):
        self.tokens = ['chr1', '2', '15', 'test', '0', '+', '4', '13', '0,128,0', '3', '4,3,3', '0,5,10']
        self.t = Transcript(self.tokens)
        self.transcript_seq = 'ATTCTGGCTA'
        self.cds_seq = 'TCTGGC'
        self.amino_acid = 'SG'
        self.chrom_seq = {'chr1': 'GTATTCTTGGACCTAA'}

    def test_sizes(self):
        """
        Make sure sizes are correct
        """
        self.assertEqual(len(self.t), len(self.transcript_seq))
        self.assertEqual(len(self.t.get_cds(self.chrom_seq)), len(self.cds_seq))
        self.assertEqual(len(self.t.get_protein_sequence(self.chrom_seq)), len(self.amino_acid))

    def test_chromosome_invalid_coordinates(self):
        """
        chromosome coordinate translations should return None if the coordinate is invalid
        in other spaces
        """
        for i in (-10, -1, 0, 1, 6, 10, 11, 16, 100):
            self.assertIsNone(self.t.chromosome_coordinate_to_mrna(i))
        for i in (-10, -1, 0, 1, 2, 3, 6, 10, 11, 15, 100):
            self.assertIsNone(self.t.chromosome_coordinate_to_cds(i))

    def test_mrna_invalid_coordinates(self):
        """
        mrna coordinate translation should return None if the coordinate is invalid
        in other spaces
        """
        for i in (-10, -1, 16, 100):
            self.assertIsNone(self.t.mrna_coordinate_to_chromosome(i))
        for i in (-10, -1, 0, 1, 9, 10, 100):
            self.assertIsNone(self.t.mrna_coordinate_to_cds(i))

    def test_cds_invalid_coordinates(self):
        """
        CDS coordinate translations should return None if the coordinate is invalid
        in other spaces
        """
        for i in (-10, -1, 6, 100):
            self.assertIsNone(self.t.cds_coordinate_to_chromosome(i))
            self.assertIsNone(self.t.cds_coordinate_to_mrna(i))

    def test_chromosome_coordinate_translations(self):
        """
        Check all possible chromosome translations for correct result
        """
        cds_result = [None, None, None, None, 0, 1, None, 2, 3, 4, None, None, 5, None, None, None]
        mrna_result = [None, None, 0, 1, 2, 3, None, 4, 5, 6, None, None, 7, 8, 9, None]
        for i in range(16):
            self.assertEqual(self.t.chromosome_coordinate_to_cds(i), cds_result[i])
            self.assertEqual(self.t.chromosome_coordinate_to_mrna(i), mrna_result[i])

    def test_mrna_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = [2, 3, 4, 5, 7, 8, 9, 12, 13, 14, None]
        cds_result = [None, None, 0, 1, 2, 3, 4, 5, None, None, None]
        for i in range(11):
            self.assertEqual(self.t.mrna_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.mrna_coordinate_to_cds(i), cds_result[i])

    def test_cds_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = [4, 5, 7, 8, 9, 12]
        mrna_result = [2, 3, 4, 5, 6, 7]
        for i in range(6):
            self.assertEqual(self.t.cds_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.cds_coordinate_to_mrna(i), mrna_result[i])

    def test_reciprocal_translations(self):
        """
        Test reciprocal translations between coordinate spaces
        """
        for i in range(16):
            tmp = self.t.chromosome_coordinate_to_mrna(i)
            # can't have reciprocal connection once None appears
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_chromosome(tmp), i)

            tmp = self.t.chromosome_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.mrna_coordinate_to_chromosome(i)
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)

            tmp = self.t.mrna_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_mrna(tmp), i)

            tmp = self.t.cds_coordinate_to_mrna(i)
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.chromosome_coordinate_to_mrna(i)
            if tmp is not None:
                tmp = self.t.mrna_coordinate_to_cds(tmp)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(i)
            if tmp is not None:
                tmp = self.t.chromosome_coordinate_to_mrna(tmp)
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(self.t.mrna_coordinate_to_cds(i))
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)

    def test_sequences(self):
        """
        Tests that the proper sequences are created from the intervals
        """
        self.assertEqual(self.t.get_mrna(self.chrom_seq), self.transcript_seq)
        self.assertEqual(self.t.get_cds(self.chrom_seq), self.cds_seq)
        self.assertEqual(self.t.get_protein_sequence(self.chrom_seq), self.amino_acid)


class NegativeStrandTranscriptTests(unittest.TestCase):
    """
    Tests the Transcript functionality of sequence_lib.
    Tests the example negative strand BED record drawn out below:
    chrom    0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
    seq      G  T  A  T  T  C  T  T  G  G  A  C  C  T  A  A
    tx       -  -  a  t  T  C  -  T  G  G  -  -  C  t  a  -
    tx.pos         9  8  7  6     5  4  3        2  1  0
    cds.pos              5  4     3  2  1        0
    """

    def setUp(self):
        self.t = Transcript(['chr1', '2', '15', 'A', '0', '-', '4', '13', '0,128,0', '3', '4,3,3', '0,5,10'])
        self.transcript_seq = 'TAGCCAGAAT'
        self.cds_seq = 'GCCAGA'
        self.amino_acid = 'AR'
        self.chrom_seq = {'chr1': 'GTATTCTTGGACCTAA'}

    def test_chromosome_invalid_coordinates(self):
        """
        chromosome coordinate translations should return None if the coordinate is invalid
        in other spaces
        """
        for i in (-10, -1, 0, 1, 6, 10, 11, 15, 16, 100):
            self.assertIsNone(self.t.chromosome_coordinate_to_mrna(i))
        for i in (-10, -1, 0, 1, 2, 3, 6, 10, 11, 13, 14, 15, 100):
            self.assertIsNone(self.t.chromosome_coordinate_to_cds(i))

    def test_mrna_invalid_coordinates(self):
        """
        mrna coordinate translation should return None if the coordinate is invalid
        in other spaces
        """
        for i in (-10, -1, 16, 100):
            self.assertIsNone(self.t.mrna_coordinate_to_chromosome(i))
        for i in (-10, -1, 0, 1, 8, 9, 10, 100):
            self.assertIsNone(self.t.mrna_coordinate_to_cds(i))

    def test_cds_invalid_coordinates(self):
        """
        CDS coordinate translations should return None if the coordinate is invalid
        in other spaces
        """
        for i in (-10, -1, 6, 100):
            self.assertIsNone(self.t.cds_coordinate_to_chromosome(i))
            self.assertIsNone(self.t.cds_coordinate_to_mrna(i))

    def test_chromosome_coordinate_translations(self):
        """
        Check all possible chromosome translations for correct result
        """
        cds_result = [None, None, None, None, 5, 4, None, 3, 2, 1, None, None, 0, None, None, None]
        mrna_result = [None, None, 9, 8, 7, 6, None, 5, 4, 3, None, None, 2, 1, 0, None]
        for i in range(16):
            self.assertEqual(self.t.chromosome_coordinate_to_cds(i), cds_result[i])
            self.assertEqual(self.t.chromosome_coordinate_to_mrna(i), mrna_result[i])

    def test_mrna_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = [14, 13, 12, 9, 8, 7, 5, 4, 3, 2, None]
        cds_result = [None, None, 0, 1, 2, 3, 4, 5, None, None, None]
        for i in range(11):
            self.assertEqual(self.t.mrna_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.mrna_coordinate_to_cds(i), cds_result[i])

    def test_cds_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = [12, 9, 8, 7, 5, 4]
        mrna_result = [2, 3, 4, 5, 6, 7]
        for i in range(6):
            self.assertEqual(self.t.cds_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.cds_coordinate_to_mrna(i), mrna_result[i])

    def test_reciprocal_translations(self):
        """
        Test reciprocal translations between coordinate spaces
        """
        for i in range(16):
            tmp = self.t.chromosome_coordinate_to_mrna(i)
            # can't have reciprocal connection once None appears
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_chromosome(tmp), i)

            tmp = self.t.chromosome_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.mrna_coordinate_to_chromosome(i)
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)

            tmp = self.t.mrna_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_mrna(tmp), i)

            tmp = self.t.cds_coordinate_to_mrna(i)
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.chromosome_coordinate_to_mrna(i)
            if tmp is not None:
                tmp = self.t.mrna_coordinate_to_cds(tmp)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(i)
            if tmp is not None:
                tmp = self.t.chromosome_coordinate_to_mrna(tmp)
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(self.t.mrna_coordinate_to_cds(i))
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)


class ComplicatedTranscript1(unittest.TestCase):
    """
    Tests the Transcript functionality part of sequence_lib.
    Tests the example complicated mrna below:
    chrom    0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20
    seq      G  T  A  T  T  C  T  T  G  G  A  C  C  T  A  A  G  C  C  T  G
    tx       -  t  a  t  -  -  t  t  G  G  -  -  -  T  A  A  -  c  c  t  -
    tx.pos      0  1  2        3  4  5  6           7  8  9     10 11 12
    cds.pos                          0  1           2  3  4
    """

    def setUp(self):
        self.tokens = ['chr1', '1', '20', 'A', '0', '+', '8', '16', '0,128,0', '4', '3,4,3,3', '0,5,12,16']
        self.t = Transcript(self.tokens)
        self.transcript_seq = 'TATTTGGTAACCT'
        self.cds_seq = 'GGTAA'
        self.amino_acid = 'G'
        self.chrom_seq = {'chr1':'GTATTCTTGGACCTAAGCCTG'}

    def test_sizes(self):
        """
        Make sure sizes are correct
        """
        self.assertEqual(len(self.t), len(self.transcript_seq))
        self.assertEqual(len(self.t.get_cds(self.chrom_seq)), len(self.cds_seq))
        self.assertEqual(len(self.t.get_protein_sequence(self.chrom_seq)), len(self.amino_acid))

    def test_chromosome_coordinate_translations(self):
        """
        Check all possible chromosome translations for correct result
        """
        cds_result = [None, None, None, None, None, None, None, None, 0, 1, None, None, None, 2, 3, 4, None, None, None, None, None]
        mrna_result = [None, 0, 1, 2, None, None, 3, 4, 5, 6, None, None, None, 7, 8, 9, None, 10, 11, 12, None]
        for i in range(21):
            self.assertEqual(self.t.chromosome_coordinate_to_cds(i), cds_result[i])
            self.assertEqual(self.t.chromosome_coordinate_to_mrna(i), mrna_result[i])

    def test_mrna_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = [1, 2, 3, 6, 7, 8, 9, 13, 14, 15, 17, 18, 19]
        cds_result = [None, None, None, None, None, 0, 1, 2, 3, 4, None, None, None]
        for i in range(13):
            self.assertEqual(self.t.mrna_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.mrna_coordinate_to_cds(i), cds_result[i])

    def test_cds_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = [8, 9, 13, 14, 15]
        mrna_result = [5, 6, 7, 8, 9]
        for i in range(5):
            self.assertEqual(self.t.cds_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.cds_coordinate_to_mrna(i), mrna_result[i])

    def test_reciprocal_translations(self):
        """
        Test reciprocal translations between coordinate spaces
        """
        for i in range(-1, 12):
            tmp = self.t.chromosome_coordinate_to_mrna(i)
            # can't have reciprocal connection once None appears
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_chromosome(tmp), i)

            tmp = self.t.chromosome_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.mrna_coordinate_to_chromosome(i)
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)

            tmp = self.t.mrna_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_mrna(tmp), i)

            tmp = self.t.cds_coordinate_to_mrna(i)
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.chromosome_coordinate_to_mrna(i)
            if tmp is not None:
                tmp = self.t.mrna_coordinate_to_cds(tmp)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(i)
            if tmp is not None:
                tmp = self.t.chromosome_coordinate_to_mrna(tmp)
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(self.t.mrna_coordinate_to_cds(i))
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)

    def test_sequences(self):
        """
        Tests that the proper sequences are created from the intervals
        """
        self.assertEqual(self.t.get_mrna(self.chrom_seq), self.transcript_seq)
        self.assertEqual(self.t.get_cds(self.chrom_seq), self.cds_seq)
        self.assertEqual(self.t.get_protein_sequence(self.chrom_seq), self.amino_acid)

    def test_get_bed(self):
        self.assertEqual(self.t.get_bed(), self.tokens)
        self.assertEqual(self.t.get_bed(new_start=1, new_stop=12),
                         ['chr1', '1', '10', 'A', '0', '+', '8', '10', '0,128,0', '2', '3,4', '0,5'])
        self.assertEqual(self.t.get_bed(new_start=19, new_stop=19),
                         ['chr1', '19', '19', 'A', '0', '+', '0', '0', '0,128,0', '1', '0', '0'])
        self.assertEqual(self.t.get_bed(new_start=1, new_stop=4),
                         ['chr1', '1', '4', 'A', '0', '+', '0', '0', '0,128,0', '1', '3', '0'])


class ComplicatedTranscript2(unittest.TestCase):
    """
    Tests the Transcript functionality part of sequence_lib.
    Tests the example negative strand complicated mrna below:
    chrom    0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20
    seq      G  T  A  T  T  C  T  T  G  G  A  C  C  T  A  A  G  C  C  T  G
    tx       -  t  a  t  -  -  t  T  G  G  -  -  -  T  A  A  -  c  c  t  -
    tx.pos      12 11 10       9  8  7  6           5  4  3     2  1  0
    cds.pos                       5  4  3           2  1  0
    """

    def setUp(self):
        self.tokens = ['chr1', '1', '20', 'A', '0', '-', '7', '16', '0,128,0', '4', '3,4,3,3', '0,5,12,16']
        self.t = Transcript(self.tokens)
        self.transcript_seq = 'AGGTTACCAAATA'
        self.cds_seq = 'TTACCA'
        self.amino_acid = 'LP'
        self.chrom_seq = {'chr1': 'GTATTCTTGGACCTAAGCCTG'}

    def test_sizes(self):
        """
        Make sure sizes are correct
        """
        self.assertEqual(len(self.t), len(self.transcript_seq))
        self.assertEqual(len(self.t.get_cds(self.chrom_seq)), len(self.cds_seq))
        self.assertEqual(len(self.t.get_protein_sequence(self.chrom_seq)), len(self.amino_acid))
        self.assertEqual(len(self.t.get_cds(self.chrom_seq)), self.t.cds_size)

    def test_chromosome_coordinate_translations(self):
        """
        Check all possible chromosome translations for correct result
        """
        cds_result = [None, None, None, None, None, None, None, 5, 4, 3, None, None, None, 2, 1, 0, None, None, None, None, None]
        mrna_result = [None, 12, 11, 10, None, None, 9, 8, 7, 6, None, None, None, 5, 4, 3, None, 2, 1, 0, None]
        for i in range(21):
            self.assertEqual(self.t.chromosome_coordinate_to_cds(i), cds_result[i])
            self.assertEqual(self.t.chromosome_coordinate_to_mrna(i), mrna_result[i])

    def test_mrna_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = [19, 18, 17, 15, 14, 13, 9, 8, 7, 6, 3, 2, 1]
        cds_result = [None, None, None, 0, 1, 2, 3, 4, 5, None, None, None, None]
        for i in range(13):
            self.assertEqual(self.t.mrna_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.mrna_coordinate_to_cds(i), cds_result[i])

    def test_cds_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = [15, 14, 13, 9, 8, 7]
        mrna_result = [3, 4, 5, 6, 7, 8]
        for i in range(5):
            self.assertEqual(self.t.cds_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.cds_coordinate_to_mrna(i), mrna_result[i])

    def test_reciprocal_translations(self):
        """
        Test reciprocal translations between coordinate spaces
        """
        for i in range(-1, 12):
            tmp = self.t.chromosome_coordinate_to_mrna(i)
            # can't have reciprocal connection once None appears
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_chromosome(tmp), i)

            tmp = self.t.chromosome_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.mrna_coordinate_to_chromosome(i)
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)

            tmp = self.t.mrna_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_mrna(tmp), i)

            tmp = self.t.cds_coordinate_to_mrna(i)
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.chromosome_coordinate_to_mrna(i)
            if tmp is not None:
                tmp = self.t.mrna_coordinate_to_cds(tmp)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(i)
            if tmp is not None:
                tmp = self.t.chromosome_coordinate_to_mrna(tmp)
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(self.t.mrna_coordinate_to_cds(i))
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)

    def test_sequences(self):
        """
        Tests that the proper sequences are created from the intervals
        """
        self.assertEqual(self.t.get_mrna(self.chrom_seq), self.transcript_seq)
        self.assertEqual(self.t.get_cds(self.chrom_seq), self.cds_seq)
        self.assertEqual(self.t.get_protein_sequence(self.chrom_seq), self.amino_acid)

    def test_get_bed(self):
        self.assertEqual(self.t.get_bed(), self.tokens)
        self.assertEqual(self.t.get_bed(new_start=4),
                         ['chr1', '6', '20', 'A', '0', '-', '7', '16', '0,128,0', '3', '4,3,3', '0,7,11'])
        self.assertEqual(self.t.get_bed(new_start=17),
                         ['chr1', '17', '20', 'A', '0', '-', '0', '0', '0,128,0', '1', '3', '0'])
        self.assertEqual(self.t.get_bed(new_start=10, new_stop=17),
                         ['chr1', '13', '16', 'A', '0', '-', '13', '16', '0,128,0', '1', '3', '0'])


class SingleExonTranscript1(unittest.TestCase):
    """
    Tests the Transcript functionality part of sequence_lib.
    Tests the example single exon mrna below:
    chrom    0  1  2  3  4  5
    seq      G  T  A  T  T  C
    tx       g  T  A  T  t  c
    tx.pos   0  1  2  3  4  5
    cds.pos     0  1  2
    """

    def setUp(self):
        self.t = Transcript(['chr1', '0', '6', 'A', '0', '+', '1', '4', '0,128,0', '1', '6', '0'])
        self.transcript_seq = 'GTATTC'
        self.cds_seq = 'TAT'
        self.amino_acid = 'Y'
        self.chrom_seq = {'chr1': 'GTATTCTTGGACCTAA'}

    def test_sizes(self):
        """
        Make sure sizes are correct
        """
        self.assertEqual(len(self.t), len(self.transcript_seq))
        self.assertEqual(len(self.t.get_cds(self.chrom_seq)), len(self.cds_seq))
        self.assertEqual(len(self.t.get_protein_sequence(self.chrom_seq)), len(self.amino_acid))

    def test_mrna_invalid_coordinates(self):
        """
        mrna coordinate translation should return None if the coordinate is invalid
        in other spaces
        """
        for i in (-10, -1, 6, 100):
            self.assertIsNone(self.t.mrna_coordinate_to_chromosome(i))
        for i in (-10, -1, 0, 4, 5, 9, 10, 100):
            self.assertIsNone(self.t.mrna_coordinate_to_cds(i))

    def test_cds_invalid_coordinates(self):
        """
        CDS coordinate translations should return None if the coordinate is invalid
        in other spaces
        """
        for i in (-10, -1, 4, 100):
            self.assertIsNone(self.t.cds_coordinate_to_chromosome(i))
            self.assertIsNone(self.t.cds_coordinate_to_mrna(i))

    def test_chromosome_coordinate_translations(self):
        """
        Check all possible chromosome translations for correct result
        """
        cds_result = [None, 0, 1, 2, None, None, None]
        mrna_result = [0, 1, 2, 3, 4, 5, None]
        for i in range(6):
            self.assertEqual(self.t.chromosome_coordinate_to_cds(i), cds_result[i])
            self.assertEqual(self.t.chromosome_coordinate_to_mrna(i), mrna_result[i])

    def test_mrna_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = [0, 1, 2, 3, 4, 5, None]
        cds_result = [None, 0, 1, 2, None, None, None]
        for i in range(6):
            self.assertEqual(self.t.mrna_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.mrna_coordinate_to_cds(i), cds_result[i])

    def test_cds_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = [1, 2, 3, None]
        mrna_result = [1, 2, 3, None]
        for i in range(4):
            self.assertEqual(self.t.cds_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.cds_coordinate_to_mrna(i), mrna_result[i])

    def test_reciprocal_translations(self):
        """
        Test reciprocal translations between coordinate spaces
        """
        for i in range(-1, 7):
            tmp = self.t.chromosome_coordinate_to_mrna(i)
            # can't have reciprocal connection once None appears
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_chromosome(tmp), i)

            tmp = self.t.chromosome_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.mrna_coordinate_to_chromosome(i)
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)

            tmp = self.t.mrna_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_mrna(tmp), i)

            tmp = self.t.cds_coordinate_to_mrna(i)
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.chromosome_coordinate_to_mrna(i)
            if tmp is not None:
                tmp = self.t.mrna_coordinate_to_cds(tmp)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(i)
            if tmp is not None:
                tmp = self.t.chromosome_coordinate_to_mrna(tmp)
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(self.t.mrna_coordinate_to_cds(i))
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)

    def test_sequences(self):
        """
        Tests that the proper sequences are created from the intervals
        """
        self.assertEqual(self.t.get_mrna(self.chrom_seq), self.transcript_seq)
        self.assertEqual(self.t.get_cds(self.chrom_seq), self.cds_seq)
        self.assertEqual(self.t.get_protein_sequence(self.chrom_seq), self.amino_acid)


class SingleExonTranscript2(unittest.TestCase):
    """
    Tests the Transcript functionality part of sequence_lib.
    Tests the example single exon mrna below:
    chrom    0  1  2  3  4  5
    seq      G  T  A  T  T  C
    tx       G  T  A  T  T  C
    tx.pos   0  1  2  3  4  5
    cds.pos  0  1  2  3  4  5
    """

    def setUp(self):
        self.t = Transcript(['chr1', '0', '6', 'A', '0', '+', '0', '6', '0,128,0', '1', '6', '0'])
        self.transcript_seq = 'GTATTC'
        self.cds_seq = self.transcript_seq
        self.amino_acid = 'VF'
        self.chrom_seq = {'chr1': 'GTATTCTTGGACCTAA'}

    def test_sizes(self):
        """
        Make sure sizes are correct
        """
        self.assertEqual(len(self.t), len(self.transcript_seq))
        self.assertEqual(len(self.t.get_cds(self.chrom_seq)), len(self.cds_seq))
        self.assertEqual(len(self.t.get_protein_sequence(self.chrom_seq)), len(self.amino_acid))

    def test_chromosome_coordinate_translations(self):
        """
        Check all possible chromosome translations for correct result
        """
        cds_result = mrna_result = [0, 1, 2, 3, 4, 5, None]
        for i in range(6):
            self.assertEqual(self.t.chromosome_coordinate_to_cds(i), cds_result[i])
            self.assertEqual(self.t.chromosome_coordinate_to_mrna(i), mrna_result[i])

    def test_mrna_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = cds_result = [0, 1, 2, 3, 4, 5, None]
        for i in range(6):
            self.assertEqual(self.t.mrna_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.mrna_coordinate_to_cds(i), cds_result[i])

    def test_cds_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = mrna_result = [0, 1, 2, 3, 4, 5, None]
        for i in range(6):
            self.assertEqual(self.t.cds_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.cds_coordinate_to_mrna(i), mrna_result[i])

    def test_reciprocal_translations(self):
        """
        Test reciprocal translations between coordinate spaces
        """
        for i in range(-1, 7):
            tmp = self.t.chromosome_coordinate_to_mrna(i)
            # can't have reciprocal connection once None appears
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_chromosome(tmp), i)

            tmp = self.t.chromosome_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.mrna_coordinate_to_chromosome(i)
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)

            tmp = self.t.mrna_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_mrna(tmp), i)

            tmp = self.t.cds_coordinate_to_mrna(i)
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.chromosome_coordinate_to_mrna(i)
            if tmp is not None:
                tmp = self.t.mrna_coordinate_to_cds(tmp)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(i)
            if tmp is not None:
                tmp = self.t.chromosome_coordinate_to_mrna(tmp)
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(self.t.mrna_coordinate_to_cds(i))
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)

    def test_sequences(self):
        """
        Tests that the proper sequences are created from the intervals
        """
        self.assertEqual(self.t.get_mrna(self.chrom_seq), self.transcript_seq)
        self.assertEqual(self.t.get_cds(self.chrom_seq), self.cds_seq)
        self.assertEqual(self.t.get_protein_sequence(self.chrom_seq), self.amino_acid)


class NoncodingTranscript(unittest.TestCase):
    """
    Tests the Transcript functionality part of sequence_lib.
    Tests the example non-coding spliced mrna below:
    chrom    0  1  2  3  4  5  6  7  8  9  10
    seq      G  T  A  T  T  C  T  T  G  G  A
    tx       g  t  a  t  -  -  t  -  g  g  a
    tx.pos   0  1  2  3        4     5  6  7
    """

    def setUp(self):
        self.t = Transcript(['chr1', '0', '11', 'A', '0', '+', '0', '0', '0,128,0', '3', '4,1,3', '0,6,8'])
        self.transcript_seq = 'GTATTGGA'
        self.cds_seq = ''
        self.amino_acid = ''
        self.chrom_seq = {'chr1': 'GTATTCTTGGACCTAA'}

    def test_sizes(self):
        """
        Make sure sizes are correct
        """
        self.assertEqual(len(self.t), len(self.transcript_seq))
        self.assertEqual(len(self.t.get_cds(self.chrom_seq)), len(self.cds_seq))
        self.assertEqual(len(self.t.get_protein_sequence(self.chrom_seq)), len(self.amino_acid))

    def test_chromosome_coordinate_translations(self):
        """
        Check all possible chromosome translations for correct result
        """
        cds_result = [None] * 12
        mrna_result = [0, 1, 2, 3, None, None, 4, None, 5, 6, 7, None]
        for i in range(12):
            self.assertEqual(self.t.chromosome_coordinate_to_cds(i), cds_result[i])
            self.assertEqual(self.t.chromosome_coordinate_to_mrna(i), mrna_result[i])

    def test_mrna_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = [0, 1, 2, 3, 6, 8, 9, 10, None]
        cds_result = [None] * 9
        for i in range(9):
            self.assertEqual(self.t.mrna_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.mrna_coordinate_to_cds(i), cds_result[i])

    def test_cds_coordinate_translations(self):
        """
        Check all possible mrna translations for correct result
        """
        chrom_result = mrna_result = [None] * 10
        for i in range(10):
            self.assertEqual(self.t.cds_coordinate_to_chromosome(i), chrom_result[i])
            self.assertEqual(self.t.cds_coordinate_to_mrna(i), mrna_result[i])

    def test_reciprocal_translations(self):
        """
        Test reciprocal translations between coordinate spaces
        """
        for i in range(-1, 12):
            tmp = self.t.chromosome_coordinate_to_mrna(i)
            # can't have reciprocal connection once None appears
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_chromosome(tmp), i)

            tmp = self.t.chromosome_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.mrna_coordinate_to_chromosome(i)
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)

            tmp = self.t.mrna_coordinate_to_cds(i)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_mrna(tmp), i)

            tmp = self.t.cds_coordinate_to_mrna(i)
            if tmp is not None:
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.chromosome_coordinate_to_mrna(i)
            if tmp is not None:
                tmp = self.t.mrna_coordinate_to_cds(tmp)
            if tmp is not None:
                self.assertEqual(self.t.cds_coordinate_to_chromosome(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(i)
            if tmp is not None:
                tmp = self.t.chromosome_coordinate_to_mrna(tmp)
                self.assertEqual(self.t.mrna_coordinate_to_cds(tmp), i)

            tmp = self.t.cds_coordinate_to_chromosome(self.t.mrna_coordinate_to_cds(i))
            if tmp is not None:
                self.assertEqual(self.t.chromosome_coordinate_to_mrna(tmp), i)

    def test_sequences(self):
        """
        Tests that the proper sequences are created from the intervals
        """
        self.assertEqual(self.t.get_mrna(self.chrom_seq), self.transcript_seq)
        self.assertEqual(self.t.get_cds(self.chrom_seq), self.cds_seq)
        self.assertEqual(self.t.get_protein_sequence(self.chrom_seq), self.amino_acid)


class PositiveStrandGenePredTranscript(PositiveStrandTranscriptTests):
    """
    Tests the Transcript functionality part of sequence_lib.
    Tests the example positive strand BED record drawn out below:
    chrom    0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
    seq      G  T  A  T  T  C  T  T  G  G  A  C  C  T  A  A
    tx       -  -  a  t  T  C  -  T  G  G  -  -  C  t  a  -
    tx.pos         0  1  2  3     4  5  6        7  8  9
    cds.pos              0  1     2  3  4        5
    """

    def setUp(self):
        self.tokens = ['A', 'chr1', '+', '2', '15', '4', '13', '3', '2,7,12', '6,10,15', '1',
                       'q2', 'cmpl', 'cmpl', '2,1,1']
        self.t = GenePredTranscript(self.tokens)
        self.transcript_seq = 'ATTCTGGCTA'
        self.cds_seq = 'TCTGGC'
        self.amino_acid = 'L'  # this transcript has a offset of 2, so the first in-frame codon is TGG
        self.chrom_seq = {'chr1': 'GTATTCTTGGACCTAA'}

    def test_sequences(self):
        """
        Tests that the proper sequences are created from the intervals
        """
        self.assertEqual(self.t.get_mrna(self.chrom_seq), self.transcript_seq)
        self.assertEqual(self.t.get_cds(self.chrom_seq), self.cds_seq)
        self.assertEqual(self.t.get_protein_sequence(self.chrom_seq), self.amino_acid)

    def test_get_gp(self):
        self.assertEqual(self.t.get_gene_pred(), self.tokens)


class NegativeStrandGenePredTranscript(NegativeStrandTranscriptTests):
    """
    Tests the Transcript functionality part of sequence_lib.
    Tests the example positive strand BED record drawn out below:
    chrom    0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
    seq      G  T  A  T  T  C  T  T  G  G  A  C  C  T  A  A
    tx       -  -  a  t  T  C  -  T  G  G  -  -  C  t  a  -
    tx.pos         0  1  2  3     4  5  6        7  8  9
    cds.pos              0  1     2  3  4        5
    """

    def setUp(self):
        self.tokens = ['A', 'chr1', '-', '2', '15', '4', '13', '3', '2,7,12', '6,10,15', '1',
                                    'q2', 'cmpl', 'cmpl', '2,2,1']
        self.t = GenePredTranscript(self.tokens)
        self.transcript_seq = 'TAGCCAGAAT'
        self.cds_seq = 'GCCAGA'
        self.amino_acid = 'Q'  # this transcript has a offset of 1, so the first in-frame codon is CAG
        self.chrom_seq = {'chr1': 'GTATTCTTGGACCTAA'}

    def test_sequences(self):
        """
        Tests that the proper sequences are created from the intervals
        """
        self.assertEqual(self.t.get_mrna(self.chrom_seq), self.transcript_seq)
        self.assertEqual(self.t.get_cds(self.chrom_seq), self.cds_seq)
        self.assertEqual(self.t.get_protein_sequence(self.chrom_seq), self.amino_acid)

    def test_get_gp(self):
        self.assertEqual(self.t.get_gene_pred(), self.tokens)


if __name__ == '__main__':
    unittest.main()
