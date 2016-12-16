#Comparative-Annotation-Toolkit

This project aims to provide a straightforward end-to-end pipeline that takes as input a HAL-format multiple whole genome alignment as well as a GFF3 file representing annotations on one high quality assembly in the HAL alignment, and produces a output GFF3 annotation on all target genomes chosen.

This pipeline is capable of running both on local cluster hardware as well as on common cloud infrastructure using the [toil](http://toil.readthedocs.io/en/latest/) workflow engine. For full runs on many genomes, a decent amount of computational effort is required. Memory usage is moderate.

#Dependencies

##Binaries
1. [Augustus](http://augustus.gobics.de/). 
2. [bedtools](http://bedtools.readthedocs.io/en/latest/).
3. [samtools](http://www.htslib.org/) (1.0 or greater).

In addition to the `augustus` binary, the following scripts from the augustus repository should also be in your path: `transMap2hints.pl`, `homGeneMapping`, `bam2hints`, `bam2wig`, `filterBam`.

##External repositories that provide binaries

1. [Kent toolkit](https://github.com/ucscGenomeBrowser/kent).
2. [HAL toolkit](https://github.com/glennhickey/hal).

Both of these repositories must be downloaded and compiled, and the resulting binaries placed on your path.

##Python libraries
Python version must be 2.7.x.

1. pyfasta.
2. toil.
3. luigi.
4. seaborn.
5. frozendict.
6. configobj.

All of the python depenencies are easily pip installable. The binary dependencies must be manually installed and placed on your path.

#Execution modes

The default mode of this pipeline will perform the following tasks:

1. Lift all annotations present in the input GFF3 to all other genomes in the alignment.
2. Filter these comparative transcripts for paralogous mappings.
3. Evaluate these transcripts for potential problems, assigning a score.
4. Produce a output annotation set as well a series of plots charting how this process went.

These steps will run reasonably fast on one machine without any need for cluster computing. However, to construct a high quality annotation set, it is recommended that the pipeline be run with as many modes of `augustus` as possible.

##AugustusTM(R)
Augustus transMap (TM) and augustus transMap-RNA-seq (TMR) are two execution modes that require significant computation but also provide a higher quality result. Based on the presence of a [hints database](#Hints database) with RNA-seq information, `CAT` will decide whether to run just `augustusTM` or `augustusTMR` in addition. `augustusTM` makes use only of the `transMap` information to try and build proper transcript models for each `transMap` transcript one at a time. This fixes small indel errors in the alignment as well as splice site shifts during evolutionary time by enforcing a full coding transcript model over the transcript. `augustusTMR` takes in additional RNA-seq information, and is especially useful at larger phylogenetic distance from the source genome.

##AugustusCGP
`augustusCGP` is the comparative mode of `augustus` recently introduced by Stefanie Köneig. This mode of `augustus` takes as input a HAL format multiple whole genome alignment and simultaneously produces *de-novo* transcript predictions in all genomes, taking into account conservation. Like `augustusTM(R)`, this process relies on a [hints database](#Hints database) that must also contain genome sequence information. The RNA-seq information is optional but will provide additional confidence on the final results. This module of the `CAT` pipeline is the only module that provides the opportunity to predict species or sub-clade specific isoforms and genes.

#Running the pipeline

This pipeline makes use of [Luigi](https://github.com/spotify/luigi) to link the various steps together. First, start the `luigid` daemon:

`luigid --background --logdir luigi_logs`

Which provides the central scheduler as well as the web UI, which can be accessed at `localhost:8082`. If you don't want to use the daemon, add the flag `--local-scheduler` to the invocation.

To run the test data, change directories to the CAT installation folder and do the following:

`export PYTHONPATH=./ && luigi --module cat RunCat --hal=test_data/vertebrates.hal --ref-genome=mm10 --config=test_config.cfg`

#Config file
A major component of producing high quality comparative annotations is making use of RNA-seq information. This information is used as hints to the `augustus` gene finding tool along with `transMap`, and is a major component of cleaning up transcript projections. This is also useful if you run the `augustusCGP` portion of the pipeline.

If you are running `augustusCGP`, then a hints database is required. At minimum, the hints database must be loaded with the genome sequence information for each target genome. Including transcript annotations for the high quality genomes in the alignment is highly recommended.

The `CAT` toolkit provides a simple `luigi` based interface for producing a hints database with the program `generate_hints_db.py`. The command line options for this program are described in the [command line options](#Command line options) section. 

To produce a hints database, you need to create a config file. The template for the config file is:

~~~~
[ANNOTATION]
Genome = /path/to/reference/gff3

[BAM]
Genome = /path/to/fofn OR /path/to/bam1.bam, /path/to/bam2.bam

[INTRONBAM]
Genome = /path/to/fofn/of/noisy/rnaseq
~~~~

Note that the BAM/INTRONBAM fields can be populated either with a comma separated list of BAMs or a single file with a line pointing to each BAM (a FOFN, or file-of-file-names). The reference sequence information will be extracted from the HAL alignment.

##RNA-seq libraries

It is **extremely** important that you use high quality RNA-seq. Libraries should be poly-A selected and paired end with a minimum read length of 75bp. If any of these are not true, it is advisable to place these libraries in the INTRONBAM field. Any genome can have a mix of BAM and INTRONBAM hints.

#Input

The pipeline has 2 required inputs if run in default mode:

1. A HAL format alignment containing your reference genome and target genomes. It is highly suggested that sufficient outgroups are incorporated. 
2. A GFF3 format annotation file for your reference genome. This annotation file should be in the format found when downloaded from [Ensembl](http://uswest.ensembl.org/info/data/ftp/index.html), or [GENCODE](http://www.gencodegenes.org/) for human/mouse. There are a few key components of this GFF3 that you need to have if you are producing your own. [See the GFF3 section](#GFF3).

If you want to run `AugustusTMR` or `AugustusCGP`, then you need to also provide a [hints database](#Hints database). 

#Modules

Each of these modules can be ran individually by changing the `luigi` invocation. 

##PrepareFiles

This module parses the GFF3 annotation input, creating a genePred format file as well as a sqlite database. In addition, sequence files for all target genomes are extracted and converted to 2bit.

This module will populate the folders `--work-dir/reference` and `--work-dir/genome_files`.

##Chaining

This step is the first precursor step to `transMap`. Pairwise genomic Kent-style chains are produced for each target genome from the designated reference. This step uses `Toil` and can be parallelized on a cluster.

This module will populate the folder `--work-dir/chaining`.

##TransMap

This step runs `transMap`. The chain files are used to project annotations present in the GFF3 from the reference genome to each target genome.

##EvaluateTransMap

This step performs the preliminary classification of `transMap` transcripts. This step populates the `TransMapEvaluation` table in the sqlite database for each target genome with the following classifiers:

1. Paralogy: The # of times this transcript was aligned.
2. AlnExtendsOffConfig: Does this alignment run off the end of a contig?
3. AlignmentPartialMap: Did this transcript not map completely?
4. AlnAbutsUnknownBases: Does this alignment have Ns immediately touching any exons?
5. AlnContainsUnknownBases: Are there any Ns within the transcript alignment?
6. Synteny: Counts the number of genes in linear order that match up to +/- 3 genes.
7. TransMapOriginalIntrons: The number of transMap introns within a wiggle distance of a intron in the parent transcript
   in transcript coordinates.
   
This module will populate the folder `--work-dir/transMap`.

##FilterTransMap

Resolves paralogs in transMap output based on MLE estimate of the distribution of alignment identities in the transMap
process. A normal distribution is fit to the -log(1 - identity) where identity != 1, which is simply a transformation of the underlying lognormal distribution. This is performed only on transcripts who map over in a 1-1 fashion. A cutoff is established as one standard deviation from the mean of this fit. 

This process is performed separately for each **transcript biotype** present in the annotation set. *This is one of many reasons why it is very important that your reference be properly biotyped!* Transcript projections whose identity are below the cutoff will be marked as `Failing`, which is a key component of the [consensus finding process](##consensus).

The output of this step populates the sqlite tables `TransMapFilterEvaluation` and `TransMapIdentityCutoffs` for each target genome. The table `TransMapIdentityCutoffs` records the identity cutoff established for each biotype.

If the `--resolve-split-genes` flag is set, this step will also resolve mappings of the same gene to multiple contigs. This process
combines information from both synteny and phylogenetic distance to determine which contig is likely the parental
contig. This should be used carefully on genomes that have a low N50. Looking at the `AlnExtendsOffConfig` field from the `EvaluateTransMap` step may help guide this decision.

This module will populate the folder `--work-dir/filtered_transMap`.

##Augustus

As [discussed above](#AugustusTM(R)), this module runs `AugustusTM(R)`. If the pipeline is ran without a hints database, only the `AugustusTM` mode will be executed. This process is one of the most computationally intensive steps, and should not be ran without a cluster.

This module will populate the folder `--work-dir/augustus`.

##AugustusCgp

As [discussed above](#augustuscgp), this module runs `AugustusCGP`. This process requires a hints database (although no RNA-seq is required, it is highly recommended). Running `AugustusCGP` is trickier than other modes. If your genomes are not closely related to an existing training set, you may need to perform logistic regression to train `AugustusCGP` before execution. A default parameter set is provided. This mode is also computationally intensive, and requires a cluster.

Each output transcript are assigned a parental gene, if possible. Parental gene assignment is done by looking to see if this transcript has at least 1 exonic base overlap with any [filtered TransMap](##filtertransmap). If the transcript overlaps more than one gene, the [Jaccard metric](http://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html) is used to try and resolve the ambiguity. If no gene stands out, this transcript is discarded.

Transcripts which are not assigned a parental gene will be considered *novel* in the [consensus finding](##consensus) step. Most often, these are the result of gene family expansion or contraction in the reference. Looking at the raw `transMap` track in the final [assembly hub](##assemblyhub) will help resolve this.

This module will populate the folder `--work-dir/augustus_cgp`.

##Hgm

`homGeneMapping` is a companion tool of `AugustusCGP`. This tool uses a HAL alignment to project RNA-seq and annotation information to target genomes. This is used to validate a splice junction in a target genome as being supported in one or more alternative genomes, as well as being supported in the reference annotation. This module populates the `*_Hgm` database table, where `*` is one of `transMap`, `augTM`, `augTMR`, or `augCGP` depending on the transcripts being evaluated. This table has three columns, `RnaSeqSupportIntronVector` and `AnnotationSupportIntronVector`. The first two columns are a result of parsing the full string, producing a integer vector of values for each splice junction representing the number of genomes in the hints database with RNA-seq support for this junction or annotation support for this junction, respectively. 

This module will populate the folder `--work-dir/hgm`.

The output of the `homGeneMapping` module has more information embedded in the output files. Each GTF format file in the above folder has a added column on the end with a string like:

`"0E-6273,1E-1524,2N:M*-1,3E-742,4E-1912,5E-1208"`

Which can be interpreted as 'species 0 had 6273 extrinsic hints (RNA-seq coverage), species 1 has 1524 extrinsic hints, species 2 (the reference) had both a non-coding (N) and coding (M) junction', and so on. The species numeric values are at the top of the file, and correlate to the species ID assigned internally in the hints database. These data can be useful if you want to dig in to a specific annotation.

##AlignTranscripts

Transcript alignment is only performed on protein coding transcripts. For `AugustusCGP`, only translation alignment is performed, using the in-frame CDS. For `transMap` and `AugustusTM(R)`
transcripts, two alignments will be performed - full mRNA and in-frame CDS. The results of these alignments are saved in the folder `--work-dir/transcript_alignment`. These alignments are used to create functional annotations of transcripts in the [EvaluateTranscripts](#evaluatetranscripts) module. 

For `AugustusCGP`, transcripts are aligned to each protein coding transcript present in their assigned parental gene.

##EvaluateTranscripts

A series of classifiers that evaluate transcript pairwise alignments for `transMap`, `AugustusTM(R)` and `AugustusCGP` output.

These classifiers are broken down into 2 groups, which will each end up as a table in the database:

\<alnMode\>\_\<txMode\>\_Metrics:

These classifiers are per-transcript evaluations based on both the transcript alignment and the genome context.

1. PercentUnknownBases: % of mRNA bases that are Ns.
2. AlnCoverage: Alignment coverage in transcript space.
3. AlnIdentity: Alignment identity in transcript space.
5. PercentMissingIntrons: Number of original introns not within a wiggle distance of any introns in the target.
6. PercentMissingExons: Do we lose any exons? Defined based on parent sequence, with wiggle room.
7. CdsStartStat: Is the CDS a complete start?
8. CdsEndStat: Is the CDS a complete stop?

\<alnMode\>\_\<txMode\>\_Evaluation:

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


Where txMode is one of transMap, augTM, augTMR, augCGP and alnMode is one of CDS or mRNA.

The evaluation tables will be loaded as tracks in the final [assembly hub](##assemblyhub).

##Consensus

The consensus finding process combines the data produced by most of the above modules to filter and produce a single gene set for each target genome.

First, each potential ortholog is assigned a score. The consensus score is a weighted average of alignment metrics and evaluations. 

~~~~
If RNA-seq is provided: 0.05 * cov + 0.5 * identity + 0.3 * rnaseq_support + 0.15 * structure score
If RNA-seq is not provided: 0.05 * cov + 0.7 * identity + 0.25 * structure score
structure score = 0.7 * percent original introns + 0.2 * percent original exons + 0.1 * evaluation score
evaluation score = evaluation score = 1 - (I(in frame stop) + I(coding indel) + I(CdsStartStat = imcpl and == cmpl in reference) + I(CdsEndStat = imcpl and == cmpl in reference)) / 4
~~~~

Each transcript is also assigned a class of `Excellent`, `Passing`, `Failing`. `Failing` is defined as transcripts whose identity is below the per-biotype cutoff established in [FilterTransMap](##filtertransmap). A transcript can be upgraded to `Excellent` by the following:

~~~~
If RNA-seq is provided: evaluation score == 1 and rnaseq_support >= 0.8
If RNA-seq is not provided: evaluation score == 1
~~~~

If all transcripts for a given source gene are `Failing`, then one longest transcript is chosen for the locus, and the gene is marked as `Failing`. Otherwise, for each transcript within the gene, the highest scoring transcript alignment is chosen. Note that this means that it is possible for the CDS-space alignments of non-CGP transcripts to be chosen. This also provides CGP the possibility to take the place of a `transMap` transcript or its derivatives, if they are more divergent. However, the algorithm gives preference first to valid mRNA space alignments. 

If `augustusCGP` was ran, transcripts which were not assigned to a parent gene are included as novel transcripts. In addition, CGP transcripts are evaluated as possible novel isoforms. Novel isoforms must have at least one splice junction that is both supported by RNA-seq and *not* supported by the reference annotation.

After consensus finding, a final output gene set is produced in both `GFF3` and `genePred` format. The `genePred` annotations also have a additional `.gp_info` file that has the additional fields described below.

The output will appear in `--output-dir/consensus`.

##Plots

A large range of plots are produced in `--output-dir/plots`. These include:

1. `cgp_novel.pdf`: The first page shows the number of genes predicted to be novel. These transcripts will have the tag `gene_biotype` tag be `'unknown_likely_coding` as we as not having a `source_gene` tag. The second page shows the number of transcripts which had novel splices and were assigned as new isoforms. These will have the `transcript_biotype` tag be `unknown_likely_coding`.
2. `completeness.pdf`: The number of genes/transcripts successfully mapped over. The top of the x axis is marked with a red line representing the amount of genes/transcripts the source annotation had.
3. `consensus_splice_support.pdf`: A paired boxplot and violinplot showing the percent of splice junctions per transcript supported by RNA-seq in any of the genomes in the database. Single exon transcripts are ignored.
4. `coverage.pdf`: A paired boxplot and violinplot that shows the overall transcript coverage in the *consensus* set. 
5. `identity.pdf`: A paired boxplot and violinplot that shows the overall transcript identity in the *consensus* set. 
6. `transmap_coverage.pdf`: A paired boxplot and violinplot that shows the overall transcript coverage in the filtered transMap output.
7.  `transmap_ identity.pdf`: A paired boxplot and violinplot that shows the overall transcript identity in the filtered transMap output.
8. `gene_failure.pdf`: The number of genes either missing or failed. Missing are genes which did not have any mappings after filtering. Failed are genes which had mappings but had no transcripts which passed filtering. One longest transcript was chosen for each of these, and they are tagged `failed_gene=True`.
9. `transcript_failure.pdf`: The number of transcripts either missing or failed. Missing are transcripts which did not have any mappings after filtering. Failed are transcripts which had an identity below the established threshold for their biotype. These were discarded unless their gene had no better candidates. These have the tag `transcript_class=Failed`. 
10. `paralogy.pdf`: Stacked bar charts of the number of alignments a given source transcript had in each target.
11. `split_genes.pdf`: The number of genes split if the `--resolve-split-genes` flag was set. The right plot shows the number of transcript alignments removed by resolution.
12. `transcript_categories.pdf`: Transcripts are categorized either excellent, passing or failing. Excellent transcripts are those well supported by RNA-seq with high identity. For protein coding transcripts, they also must have no frame-shifting indels and valid CDS ends. Passing are transcripts less supported by RNA-seq with high identity. Failing transcripts have identity below the established threshold for their biotype.
13. `transcript_modes.pdf`: The number of modes that supported a given transcript projection. Applies only to protein coding transcripts. Combinations below 1% of the total are grouped into other. Shows how many modes of analysis agree on a transcript, as well as the rate that augustusTMR rescues missing information.


###GFF3 tags:

1. `source_transcript`: The name of the parent transcript, if it exists.
2. `source_gene`: The name of the parent gene, if it exists.
3. `source_gene_common_name`: The common name of the parent gene, if it is different from the source gene.
4. `transcript_mode`: The name of the mode of operation that generated this transcript.
5. `alternative_source_transcripts`: A comma separated list of alternate IDs for this transcript.
6. `failed_gene`: This transcript is the single representative for a failed transcript.
7. `transcript_class`: One of failed, passing, excellent, novel
8. `paralogy`: The number of paralogs that were mapped over when transMap mapped this transcript
9. `gene_biotype`: gene biotype. If this is a novel gene, this will be `unknown_likely_coding`.
10. `transcript_biotype`: transcript biotype. If this is a novel gene, this will be `unknown_likely_coding`.
11. `paralog_status`: `confident` if the paralog was confidently resolved, `not_confident` if it was not.
12. `alternative_source_transcripts`: Other possible transcripts, if this was collapsed as the result of deduplication
13. `gene_alternate_contigs`: If the `--resolve-split-genes` flag was set, contigs that this gene was also found on are
    comma separated in this tag.
    
For `GFF3` output, the consensus score is in the score field. For `.gp_info`, it is a column.


#Command line options

## CAT / cat.py

`--hal`: Input HAL alignment file.

`--ref-genome`: Reference genome sequence name. Must be present in HAL.

`--out-dir`: Output directory. Defaults to `./cat_output`.

`--work-dir`: Working directory. Defaults to `./cat_work`. Stores all the intermediate files as well as the `toil` jobStore. Can be removed after completion (but not if you want to re-do any steps).

`--target-genomes`: List of genomes to use. If not set, all non-reference genomes in the HAL are used.

`--workers`: Number of local cores to use. If running `toil` in singleMachine mode, care must be taken with this value.

`--augustus`: Run AugustusTM(R)? Running AugustusTMR is determined if a hints database is provided that contains RNA-seq information.

`--augustus-species`: What Augustus species do we want to use? See the Augustus manual for more information. For mammals, human is a good choice, and this is the default value.

`--augustus-cgp`: Run AugustusCGP?

`--cgp-param`: Parameters file after training CGP on the alignment. Defaults to the default parameters.

`--resolve-split-genes`: Run split gene resolution? Not a good idea if N50 is low.

`--cgp-splice-support`: Percent of splice junctions in a CGP prediction that must be supported by RNA-seq in order for the transcript to be included. Defaults to 0.8.

`--cgp-num-exons`: Number of exons a CGP prediction must have before included. This value is particularly important with noisy RNA-seq.

See below for `toil` options shared with the hints database pipeline.


##BuildHints / generate\_hints\_db.py

Constructing the hints database is its own Luigi pipeline, which uses a separate config file. See the section near the top for a description of the config file. 

`--config`: Path to config file.

`--hal`: Path to HAL file.

`--augustus-hints-db`: Path to output sqlite database that will be passed to `CAT`.

`--work-dir`: Defaults to `./hints_work`. Stores the genomes and the hints GFFs. Can be cleaned up after.

`--workers`: Number of local cores to use. If running `toil` in singleMachine mode, care must be taken with this value.

##Toil

The remaining options are passed directly along to `toil`:

`--batchSystem`: Batch system to use. Defaults to singleMachine. If running in singleMachine mode, no cluster jobs will be submitted. In addition, care must be taken to balance the `--maxCores` field with the `--workers` field with the toil resources in `luigi.cfg`. Basically, you want to make sure that your # of toil resources multiplied by your `--maxCores` is fewer than the total number of system cores you want to use. However, I **highly** recommend using a non-local batch system. See the toil documentation for more.

`--maxCores`: The number of cores each `toil` module will use. If submitting to a batch system, this limits the number of concurrent submissions.


