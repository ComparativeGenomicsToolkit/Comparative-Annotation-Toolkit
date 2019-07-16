"""
 file:    hgm.py
 descr.:  runs homGeneMapping from the Augustus package
          for cross-species comparison of gene sets. HomGeneMapping uses
          the intron hints from the database to retrieve on a per-transcript
          basis information, which introns have RNA-Seq splice junctions (SJ) supported in which of the
          input genomes. It essentially adds the "hgm_info" string to the last column of the gtf, 

          hgm_info "16E-27,0E-13,1,2E-13,3E-1,4,5,6E-48,7E-30,8E-1,9,10E-1,11,12E-19,13E-27,14E-46,15E-6,17E-4";
          
          that encodes genome name, type of evidence and multiplicity, e.g.
          in the example above, the intron has RNA-Seq SJ support in species 16 (with mult=26),
          in species 0 (with mult=13), in species 2 (with mult=13), in species 3 (with mult=1), etc.
          The header in the gtf, gives a list of species numbers and corresponding names, e.g.
          
          # 0     129S1_SvImJ
          # 1     AKR_J
          # 2     A_J
          ...
          # 17    WSB_EiJ

 authors: Stefanie Koenig, Ian Fiddes
"""
import os
import collections
import pandas as pd
import tools.fileOps
import tools.misc
import tools.procOps
import tools.nameConversions
import tools.transcripts
import tools.hintsDatabaseInterface


def hgm(args):
    """
    Main entry function for hgm module. Runs homGeneMapping after extracting exons and non-coding introns from the
    annotation. Exons and non-coding introns are not stored in the database to make predictions go faster.
    When we are running homGeneMapping on non-CGP data, we have no new annotations for the reference, so we pass a
    dummy GTF.

    :param args: dictionary of arguments from CAT
    :return: a dictionary with one gtf file per genome
    """
    tools.fileOps.ensure_dir(args.gtf_out_dir)
    # keep track of the GTFs we are generating to remove later
    supplementary_gffs = []

    with tools.fileOps.TemporaryFilePath() as gtf_fofn, tools.fileOps.TemporaryDirectoryPath() as temp_dir:
        with open(gtf_fofn, 'w') as outf:
            for genome, gtf in args.in_gtf.iteritems():
                if genome != args.ref_genome:
                    supplementary_gff = create_supplementary_gff(args.hints_db, gtf, genome)
                else:
                    supplementary_gff = create_supplementary_gff(args.hints_db, gtf, genome, args.annotation_gp)
                if os.environ.get('CAT_BINARY_MODE') == 'singularity':
                    tools.fileOps.print_row(outf, [genome] + map(tools.procOps.singularify_arg, [gtf, supplementary_gff]))
                else:
                    tools.fileOps.print_row(outf, [genome, gtf, supplementary_gff])
                supplementary_gffs.append(supplementary_gff)
            if args.ref_genome not in args.in_gtf:  # we are not running CGP, and so have no GTF for the reference
                dummy_gtf = tools.fileOps.get_tmp_file()
                tools.fileOps.touch(dummy_gtf)
                supplementary_gff = create_supplementary_gff(args.hints_db, args.annotation_gtf, args.ref_genome,
                                                             args.annotation_gp)
                if os.environ.get('CAT_BINARY_MODE') == 'singularity':
                    tools.fileOps.print_row(outf, [args.ref_genome] + map(tools.procOps.singularify_arg, [dummy_gtf, supplementary_gff]))
                else:
                    tools.fileOps.print_row(outf, [args.ref_genome, dummy_gtf, supplementary_gff])
                supplementary_gffs.append(supplementary_gff)
            else:
                dummy_gtf = None

        cmd = ['homGeneMapping',
               '--halfile={}'.format(args.hal),
               '--dbaccess={}'.format(args.hints_db),
               '--gtfs={}'.format(gtf_fofn),
               '--outdir={}'.format(args.gtf_out_dir),
               '--tmpdir={}'.format(temp_dir),
               '--cpu={}'.format(args.hgm_cpu)]
        tools.procOps.run_proc(cmd, stdout='/dev/null')

    # cleanup
    for gff in supplementary_gffs:
        os.remove(gff)
    if dummy_gtf is not None:
        os.remove(dummy_gtf)


def create_supplementary_gff(hints_db, in_gtf, genome, annotation_gp=None):
    """
    Creates the supplementary GFF which contains exon hints derived from the database as well as non-coding introns
    and all exons if annotation_gp is passed.
    :param hints_db: path to the hints database
    :param in_gtf: GTF file for this genome. If we are not doing this on CGP results, and this is the reference genome,
    this will be the annotation GTF.
    :param genome: current genome
    :param annotation_gp: annotation genePred, if we have one
    :return: file path
    """
    hints = extract_exon_hints(hints_db, in_gtf, genome)
    if annotation_gp is not None:
        hints.extend(extract_exons_non_coding_introns(annotation_gp))
    tmp_path = tools.fileOps.get_tmp_file()
    with open(tmp_path, 'w') as outf:
        tools.fileOps.print_rows(outf, hints)
    # sort and merge hints on the same intervals
    cmd = [['sort', '-n', '-k4,4', tmp_path],
           ['sort', '-s', '-n', '-k5,5'],
           ['sort', '-s', '-k3,3'],
           ['sort', '-s', '-k1,1'],
           ['join_mult_hints.pl']]
    supplementary_gff_path = tools.fileOps.get_tmp_file(suffix='gff')
    tools.procOps.run_proc(cmd, stdout=supplementary_gff_path)
    os.remove(tmp_path)
    return supplementary_gff_path


def extract_exons_non_coding_introns(annotation_gp):
    """
    Extracts a GTF file with the exons and the non-coding splice junctions for use by homGeneMapping to recognize if
    a exon or a non-coding intron junction is supported by RNA-seq.

    :param annotation_gp: genePred of annotation
    :return: list of gff-formatted lists
    """
    hints = []
    for tx in tools.transcripts.gene_pred_iterator(annotation_gp):
        for intron in tx.intron_intervals:
            if not intron.subset(tx.coding_interval):
                r = [tx.chromosome, 'tmp', 'intron', intron.start + 1, intron.stop, '.', tx.strand, '.', 'source=N']
                hints.append(r)
        for exon in tx.exon_intervals:
            r = [tx.chromosome, 'tmp', 'exon', exon.start + 1, exon.stop, '.', tx.strand, '.', 'source=M']
            hints.append(r)
    return hints


def extract_exon_hints(hints_db, in_gtf, genome):
    """
    The hints database only contains exonpart hints. For homGeneMapping, we want to merge these intervals into
    exon hints, taking the weighted average of the scores as the overall expression. To do this, I make use of bedtools.

    After I extract the merged averaged exon intervals, I intersect this with the given annotation set for this
    comparison. I then provide these intervals as exon-hints after merging the hints.

    :param hints_db: Path to the hints database
    :param in_gtf: GTF file for this genome
    :param genome: Genome in question
    :return: list of gff-formatted lists
    """
    # extract all exonpart hints from the database
    speciesnames, seqnames, hints, featuretypes, session = tools.hintsDatabaseInterface.reflect_hints_db(hints_db)
    hints_file = tools.fileOps.get_tmp_file()
    # there is a weird bug in sqlalchemy where sometimes empty genomes hang indefinitely.
    # Either way, we don't want to attempt to extract exon hints from genomes without exon (wiggle) hints anyways
    if tools.hintsDatabaseInterface.genome_has_no_wiggle_hints(hints_db, genome):
        return []
    with open(hints_file, 'w') as outf_h:
        wiggle_iter = tools.hintsDatabaseInterface.get_wiggle_hints(genome, speciesnames, seqnames, hints, session)
        for seqname, start, end, score in wiggle_iter:
            outf_h.write('\t'.join(map(str, [seqname, start, end, score])) + '\n')
    # merge exonpart hints, averaging the coverage
    merged_hints_file = tools.fileOps.get_tmp_file()
    cmd = ['bedtools', 'merge', '-i', hints_file, '-c', '4', '-o', 'mean']
    tools.procOps.run_proc(cmd, stdout=merged_hints_file, stderr='/dev/null')
    # overlap the merged exons with the given GTF, producing a final set.
    tmp_bed = tools.fileOps.get_tmp_file()
    cmd = [['grep', '-P', '(\texon\t|\tCDS\t)', in_gtf],  # exons or CDS only
           ['cut', '-d', '\t', '-f', '1,4,5']]  # slice into BED-like format with GTF intervals
    tools.procOps.run_proc(cmd, stdout=tmp_bed)
    # sort the BED
    tools.procOps.run_proc(['bedSort', tmp_bed, tmp_bed])
    # intersect with hints and retain scores
    cmd = [['bedtools', 'intersect', '-a', tmp_bed, '-b', merged_hints_file, '-f', '0.8', '-wa', '-wb'],
           # bedtools reports both entire A and entire B if at least 80% of A overlaps a B
           ['cut', '-d', '\t', '-f', '1,2,3,7']]  # retain the A positions with the B score
    # these BED-like records are actually GFF intervals with 1-based starts and closed intervals
    bed_plus_1 = tools.procOps.call_proc_lines(cmd)

    hints = []
    for line in bed_plus_1:
        chrom, start, end, score = line.split()
        tags = 'pri=3;source=E;mult={}'.format(int(round(float(score))))
        hints.append([chrom, 'tmp', 'exon', start, end, '.', '.', '.', tags])
    os.remove(hints_file)
    os.remove(merged_hints_file)
    return hints


def parse_hgm_gtf(hgm_out, genome):
    """
    parses the hgm output gtfs and creates for each transcript a string with the intron support counts
    For now, we just count for each of the introns in the transcript, the number of species
    in which it has RNA-Seq SJ support (number of "E" in the 'hgm_info' string).
    But this can be changed later on, e.g. using also the multiplicities, or the presence of the intron
    in (one of) the reference annotation(s) (if "M" is in the 'hgm_info' string, then it is an annotated intron)

    We calculate this as both the in-species and all-species vectors.
    """
    def calculate_annot_support(intron_info, cds_info, exon_info):
        intron_annot = ','.join(map(str, [x.count('M') + x.count('N') for x in intron_info]))
        cds_annot = ','.join(map(str, [x.count('M') for x in cds_info]))
        exon_annot = ','.join(map(str, [x.count('M') for x in exon_info]))
        assert len(intron_annot) + 2 == len(exon_annot) or len(intron_annot) == 0 and len(exon_annot) == 1, \
            (len(intron_annot), len(exon_annot), aln_id)
        return [intron_annot, cds_annot, exon_annot]

    def calculate_all_species(intron_info, exon_info):
        intron_rna = ','.join(map(str, [x.count('E') + x.count('PB') for x in intron_info]))
        exon_rna = ','.join(map(str, [x.count('E') + x.count('PB') for x in exon_info]))
        return [intron_rna, exon_rna]

    def calculate_in_species(intron_info, exon_info, species_id):
        def parse_entry(entry, species_id):
            recs = entry.split(',')
            for x in recs:
                if x.startswith(species_id):
                    return x[len(species_id):]
            return ''

        intron_rna = ','.join(map(str, [parse_entry(x, species_id).count('E') +
                                        parse_entry(x, species_id).count('PB') for x in intron_info]))
        exon_rna = ','.join(map(str, [parse_entry(x, species_id).count('E') +
                                      parse_entry(x, species_id).count('PB') for x in exon_info]))
        return [intron_rna, exon_rna]

    intron_lines = []
    cds_lines = []
    exon_lines = []
    species_map = {}
    seen_lines = set()  # filter out duplicate lines. Happens in CGP/PB.
    with open(hgm_out) as infile:
        for line in infile:
            if line in seen_lines:
                continue
            seen_lines.add(line)
            if line.startswith('#') and line != '###\n':
                _, species_id, species = line.split()
                species_map[species] = species_id
            if '\tintron\t' in line:
                intron_lines.append(line.rstrip().split('\t')[-1])
            elif '\tCDS\t' in line:
                cds_lines.append(line.rstrip().split('\t')[-1])
            elif '\texon\t' in line:
                exon_lines.append(line.rstrip().split('\t')[-1])

    species_id = species_map[genome]
    # make use of the sorted nature of the input GTFs to create a ordered vector
    d = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))
    for mode, group in zip(*[['intron', 'cds', 'exon'], [intron_lines, cds_lines, exon_lines]]):
        for attr_line in group:
            attributes = tools.misc.parse_gtf_attr_line(attr_line)
            d[attributes['gene_id']][attributes['transcript_id']][mode].append(attributes['hgm_info'])

    # convert to dataframe, switching the list to a comma separated string
    dd = []
    for gene_id in d:
        for aln_id in d[gene_id]:
            intron_info = d[gene_id][aln_id]['intron']
            cds_info = d[gene_id][aln_id]['cds']
            exon_info = d[gene_id][aln_id]['exon']
            if tools.nameConversions.aln_id_is_denovo(aln_id):
                tx_id = aln_id
            else:
                tx_id = tools.nameConversions.strip_alignment_numbers(aln_id)
            all_species_vectors = calculate_all_species(intron_info, exon_info)
            in_species_vectors = calculate_in_species(intron_info, exon_info, species_id)
            annot_support_vectors = calculate_annot_support(intron_info, cds_info, exon_info)
            dd.append([gene_id, tx_id, aln_id] + all_species_vectors + in_species_vectors + annot_support_vectors)

    df = pd.DataFrame(dd)
    df.columns = ['GeneId', 'TranscriptId', 'AlignmentId',
                  'AllSpeciesIntronRnaSupport', 'AllSpeciesExonRnaSupport',
                  'IntronRnaSupport', 'ExonRnaSupport',
                  'IntronAnnotSupport', 'CdsAnnotSupport', 'ExonAnnotSupport']
    df = df.set_index(['GeneId', 'TranscriptId', 'AlignmentId'])
    return df
