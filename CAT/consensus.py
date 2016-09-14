"""
Generates consensus gene set.

This module takes as input the genePreds produced by transMap, AugustusTM(R) and AugustusCGP and generates a consensus
of these, producing a filtered gene set.

This process relies on a combination of metrics and evaluations loaded to a sqlite database by the classify module.

Transcript scoring functions:

structure score =   0.5 * (1 - # missing introns / # parent introns)
                  + 0.45 * (1 - # missing exons / # parent exons)
                  + 0.05 * # exon gain / # parent exons
This function is a weighted average of the structural changes seen in the alignments. Missing original introns is
weighted slightly higher than losing exons which is slightly higher than gaining exons because gaining exons may be
a real change. We multiply by 100 to put this score on the same scale as the base transcript score.
# parent introns/exons will be adjusted for CDS alignments to only include coding introns/exons.
Exists on the range(0, 100)

evaluation penalties = 2 * I(in frame stop not in ref)
                         + I(coding indel)
                         + I(CdsStartStat = imcpl and == cmpl in reference)
                         + I(CdsEndStat = imcpl and == cmpl in reference)
This function uses classifiers in a binary fashion as indicator functions to penalize the final transcript score for
problems. Only applies to coding transcripts. This function does not penalize transcripts who are messed up in the
reference, as this would bias towards incorrect Augustus corrections.
Exists in the range (0, 5) where 0 means no penalties and 5 is the maximum penalty.

evaluation bonus = 2 * I(CdsStartStat == cmpl) + 2 * I(CdsEndStat == cmpl)
This function looks at the CdsStartStat/CdsEndStat and provides up to 2 bonus points for complete CDS ends.
Exists in [0, 2, 4]


consensus score = 50 * (1 - badness) + 25 * structure score + 25 * # supported junctions
Exists in the range (0, 100). Higher is better.

coding consensus score = consensus score + cds bonus + evaluation bonus - evaluation penalties
For coding transcripts consensus.py modifies the consensus score by removing as much as 5 and adding as much as 4

Choosing CDS vs mRNA:

align_transcripts.py aligns transMap/AugustusTM(R) transcripts both in both CDS space and mRNA space.
For CGP transcripts, only CDS space alignments are performed. Thus, we do not have a fair
comparison if we use the mRNA alignments. Since CDS sequences are more conserved, it can be expected that if alignments
are evaluated all together that the CDS alignments will have the highest scores unless something has gone seriously
wrong (often, incomplete mapping particularly for genes with long UTRs). In those cases, the mRNA alignment will win.
For this reason, all scores are considered simultaneously.

The filter requirements below ensure that the results are at least of sufficient quality.

Filtering transcripts:

Before scoring, the transcripts are filtered based on a series of minimum cutoffs. These are:
AlignmentIdentity: mRNA: >70% CDS: >80%
AlignmentCoverage: mRNA: >50% CDS: >90%
Intron Inequality: 2 * # missing original introns <= # reference introns - 1 or # reference introns < 5
PercentUnknownBases: mRNA: <5% CDS: <1%

These filters throw out low quality alignments via the identity/coverage filter as well as the percent unknown bases.
The intron inequality is an important filter to preventing retroposed pseudogenes from being assigned as the parent.

If the filter process throws out all transcripts for a gene, then one transcript with the lowest badness score
will be kept to represent the locus.

Paralogy: If a transMap transcript mapped to more than one place, then we make use of the synteny score to decide
which locus is more correct. Generally, paralogous alignments are one of two things: actual paralogy, and alignment
chain breaking rearrangements that lead to fractured transcripts. In the first case, the synteny score should help. In
the second case, the synteny score will be generally the same. If this happens, then hopefully AUGUSTUS in one of its
forms rescued the transcript.

"""
