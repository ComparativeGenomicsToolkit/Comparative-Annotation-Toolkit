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
import tools.procOps
import tools.nameConversions
import tools.transcripts


def hgm(args):
    """
    Main entry function for hgm toil pipeline
    :param args: dictionary of arguments from CAT
    :return: a dictionary with one gtf file per genome
    """
    tools.fileOps.ensure_dir(args.gtf_out_dir)
    with tools.fileOps.TemporaryFilePath() as gtf_fofn, tools.fileOps.TemporaryDirectoryPath() as temp_dir:
        with open(gtf_fofn, 'w') as outf:
            tools.fileOps.print_rows(outf, args.in_gtf.items())
            if args.ref_genome not in args.genomes:  # create a dummy GTF for the reference
                fake_gtf = tools.fileOps.get_tmp_file()
                tools.fileOps.print_row(outf, [args.ref_genome, fake_gtf])
            else:
                fake_gtf = None

        non_coding_gtf = extract_non_coding_introns(args.annotation_gp)

        cmd = ['homGeneMapping',
               '--halfile={}'.format(args.hal),
               '--dbaccess={}'.format(args.hints_db),
               '--gtfs={}'.format(gtf_fofn),
               '--outdir={}'.format(args.gtf_out_dir),
               '--tmpdir={}'.format(temp_dir),
               '--cpu={}'.format(args.num_cpu)]
        tools.procOps.run_proc(cmd, stdout='/dev/null')

        if fake_gtf is not None:
            os.remove(fake_gtf)
        os.remove(non_coding_gtf)


def extract_non_coding_introns(annotation_gp):
    """
    Extracts a GTF file with only the non-coding splice junctions for use by homGeneMapping to recognize if
    a intron junction is supported by RNA-seq.

    This is done through the complement of intersection. First, we extract the CDS BED records for all coding
    transcripts. Then, we convert the reference genePred to GTF. We then intersect the first BED with this GTF with the
    -v flag set, returning only non-coding intervals. Finally, we filter these for intron records and add src=N to
    the last column.

    :param annotation_gp: genePred of annotation
    :return: file path
    """
    out_gtf = tools.fileOps.get_tmp_file()
    with open(out_gtf, 'w') as outf:
        for tx in tools.transcripts.gene_pred_iterator(annotation_gp):
            for intron in tx.intron_intervals:
                if not intron.subset(tx.coding_interval):
                    r = [tx.chromosome, 'tmp', 'intron', intron.start + 1, intron.stop, '.', tx.strand, '.', 'source=N']
                    tools.fileOps.print_row(outf, r)
    return out_gtf


def parse_gtf_attr_line(attr_line):
    """parse a GTF attributes line"""
    attr_line = [x.split(' ') for x in attr_line.replace('"', '').split('; ')]
    attr_line[-1][-1] = attr_line[-1][-1].replace(';', '')
    attrs = dict(attr_line)
    return attrs


def parse_hgm_gtf(hgm_out):
    """
    parses the hgm output gtfs and creates for each transcript a string with the intron support counts
    For now, we just count for each of the introns in the transcript, the number of species
    in which it has RNA-Seq SJ support (number of "E" in the 'hgm_info' string).
    But this can be changed later on, e.g. using also the multiplicities, or the presence of the intron
    in (one of) the reference annotation(s) (if "M" is in the 'hgm_info' string, then it is an annotated intron)
    """
    d = collections.defaultdict(lambda: collections.defaultdict(list))

    with open(hgm_out) as infile:
        # get last column of all intron lines
        intron_lines = [i.strip().split('\t')[-1] for i in infile if "\tintron\t" in i]
        for attr_line in intron_lines:
            attributes = parse_gtf_attr_line(attr_line)
            d[attributes['gene_id']][attributes['transcript_id']].append(attributes['hgm_info'])

    # convert to dataframe, switching the list to a comma separated string
    dd = []
    for gene_id in d:
        for aln_id, hgm_info in d[gene_id].iteritems():
            tx_id = tools.nameConversions.strip_alignment_numbers(aln_id)
            rnaseq_vector = ','.join(map(str, [x.count('E') for x in hgm_info]))
            annotation_vector = ','.join(map(str, [x.count('M') for x in hgm_info]))
            dd.append([gene_id, tx_id, aln_id, rnaseq_vector, annotation_vector])

    df = pd.DataFrame(dd)
    df.columns = ['GeneId', 'TranscriptId', 'AlignmentId', 'RnaSeqSupportIntronVector', 'AnnotationSupportIntronVector']
    df = df.set_index(['GeneId', 'TranscriptId', 'AlignmentId'])
    return df
