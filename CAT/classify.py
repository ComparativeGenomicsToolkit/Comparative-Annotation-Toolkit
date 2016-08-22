"""
A series of classifiers that evaluate transMap, AugustusTMR and AugustusCGP output.

These classifiers are broken down into 3 groups, which will each end up as a table in the database:

Alignment:

These classifiers apply only to the transMap alignments, and measure how well we mapped over this region:
1. Paralogy: The # of times this transcript was aligned
2. AlnExtendsOffConfig: Does this alignment run off the end of a contig?
3. AlignmentPartialMap: Did this transcript not map completely?
4. AlnAbutsUnknownBases: Does this alignment have Ns immediately touching any exons?
5. AlnContainsUnknownBases: Are there any Ns within the transcript alignment?
6. LongAlignment: Did this transcript align in a insanely long fashion? Indicative of paralogy problems.

Metrics:

These classifiers are per-transcript evaluations based on both the transcript alignment and the genome context.
1. PercentUnknownBases: % of mRNA bases that are Ns.
2. AlignmentCoverage: Alignment coverage in transcript space.
3. AlignmentIdentity: Alignment identity in transcript space.
4. NumberMissingOriginalIntrons: Number of original introns not within a wiggle distance of any introns in the target.
5. Goodness: A measure of how good the alignment is related to Jim Kent's badness score.

Evaluation:

These classifiers are per-transcript evaluations based on the transcript alignment.
Unlike the other two tables, this table stores the actual location of the problems (in genome coordinates) as a
BED-like format. In cases where there are multiple problems, they will be additional rows.
1. CodingInsertion: Do we have any frame-shifting coding insertions?
2. CodingDeletion: Do we have any frame-shifting coding deletions?
3. CodingMult3Insertion: Do we have any mod3 coding insertions?
4. CodingMult3Deletion: Do we have any mod3 coding deletions?
6. ExonGain: Do we gain any exons? Defined as having a continuous block of sequence with no alignment that is spliced.
7. ExonLoss: Do we lose any exons? Defined based on parent sequence, with wiggle room.
8. InFrameStop: Are there any in-frame stop codons?
9. CdsUnknownSplice: Do any CDS introns have an unknown splicing pattern?
10. UtrUnknownSplice: Do any UTR introns have an unknown splicing pattern?
11. NotOriginalCdsStart: Are the thickStarts the same?
12. NotOriginalCdsStop: Are the thickStops the same?
13. NotOriginalStart: Does transcription start at the same place?
14. NotOriginalStop: Does transcription end at the same place?

The Metrics and Evaluation groups will have multiple tables for each of the input methods used:
1) transMap
2) augTm
3) augTmr
4) augCgp

"""
import itertools
import math
import pandas as pd
import tools.transcripts
import tools.psl
import tools.fileOps
import tools.dataOps
import tools.bio
import tools.mathOps
import tools.fastaAlignment

# hard coded variables
long_transcript = 3 * 10 ** 6


