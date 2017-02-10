"""
Generate a hints database file from RNAseq/IsoSeq alignments for AugustusTMR/AugustusCGP.
"""
import collections
import itertools
import os
import logging

import pyfasta
import pysam
from toil.fileStore import FileID
from toil.common import Toil
from toil.job import Job

import tools.dataOps
import tools.fileOps
import tools.mathOps
import tools.misc
import tools.procOps
import tools.toilInterface
import tools.transcripts
import tools.hal
from exceptions import UserException

logger = logging.getLogger(__name__)


def hints_db(hints_args, toil_options):
    """
    Entry point for hints database Toil pipeline.
    """
    def validate_import_bam(t, bam_path, fasta_sequences, genome):
        validate_bam_fasta_pairs(bam_path, fasta_sequences, genome)
        return [FileID.forPath(t.importFile('file://' + bam_path), bam_path),
                FileID.forPath(t.importFile('file://' + bam_path + '.bai'), bam_path + '.bai')]

    fasta = pyfasta.Fasta(hints_args.fasta)
    fasta_sequences = {(x.split()[0], len(fasta[x])) for x in fasta.keys()}
    with Toil(toil_options) as t:
        if not t.options.restart:
            # load the RNA-seq data, if we have any
            bam_file_ids = {'BAM': {}, 'INTRONBAM': {}}
            for dtype in ['BAM', 'INTRONBAM']:
                if hints_args.genome not in hints_args.cfg[dtype]:
                    continue
                for bam_path in hints_args.cfg[dtype][hints_args.genome]:
                    bam_file_ids[dtype][os.path.basename(bam_path)] = validate_import_bam(t, bam_path,
                                                                                          fasta_sequences,
                                                                                          hints_args.genome)
                    logger.info('{} {} ({}) is valid.'.format(dtype, os.path.basename(bam_path), hints_args.genome))

            # load the IsoSeq data, if we have any
            iso_seq_file_ids = []
            if hints_args.genome in hints_args.cfg['ISO_SEQ_BAM']:
                for bam_path in hints_args.cfg['ISO_SEQ_BAM'][hints_args.genome]:
                    validate_bam_fasta_pairs(bam_path, fasta_sequences, hints_args.genome)
                    iso_seq_file_ids.append(validate_import_bam(t, bam_path, fasta_sequences, hints_args.genome))

            if hints_args.annotation is None:
                annotation_file_id = None
            else:
                annotation_file_id = FileID.forPath(t.importFile('file://' + hints_args.annotation),
                                                    hints_args.annotation)
            input_file_ids = {'bams': bam_file_ids, 'iso_seq_bams': iso_seq_file_ids, 'annotation': annotation_file_id}
            logger.info('{} has {} valid intron-only BAMs, {} valid BAMs and {} valid IsoSeq BAMs. '
                        'Beginning Toil hints pipeline.'.format(hints_args.genome, len(bam_file_ids['INTRONBAM']),
                                                                len(bam_file_ids['BAM']),
                                                                len(iso_seq_file_ids)))
            disk_usage = tools.toilInterface.find_total_disk_usage(input_file_ids)
            job = Job.wrapJobFn(setup_hints, input_file_ids, iso_seq_file_ids, disk=disk_usage)
            combined_hints = t.start(job)
        else:
            logger.info('Restarting Toil hints pipeline for {}.'.format(hints_args.genome))
            combined_hints = t.restart()
        tools.fileOps.ensure_file_dir(hints_args.hints_path)
        t.exportFile(combined_hints, 'file://' + hints_args.hints_path)


def setup_hints(job, input_file_ids, iso_seq_file_ids):
    """
    Generates hints for a given genome with a list of BAMs. Will add annotation if it exists.

    Pipeline structure:
      filter_bam    cat_sort_bams    (build_intron_hints, build_exon_hints)
        ^  V           ^ V                 ^  V
    setup_hints -> merge_bams -------> build_hints -------> cat_hints
          V generate_annotation_hints ------------------------^
          V generate_iso_seq_hints ---------------------------^
    Each main step (filter_bam, cat_sort_bams, build_intron_hints, build_exon_hints) are done on a subset of references
    that are then combined at the cat_hints step.
    """
    # RNA-seq hints
    filtered_bam_file_ids = {'BAM': collections.defaultdict(list), 'INTRONBAM': collections.defaultdict(list)}
    for dtype, bam_dict in input_file_ids['bams'].iteritems():
        if len(bam_dict) == 0:
            continue
        # Since BAMs are valid, we can assume that they all share the same header
        bam_file_id, bai_file_id = bam_dict.values()[0]
        bam_path = job.fileStore.readGlobalFile(bam_file_id)
        sam_handle = pysam.Samfile(bam_path)
        # triple disk usage to deal with name sorted bam
        disk_usage = tools.toilInterface.find_total_disk_usage([bam_file_id, bai_file_id]) * 3
        # generate reference grouping that will be used downstream until final cat step
        grouped_references = [tuple(x) for x in group_references(sam_handle)]
        for original_path, (bam_file_id, bai_file_id) in bam_dict.iteritems():
            for reference_subset in grouped_references:
                j = job.addChildJobFn(namesort_bam, bam_file_id, bai_file_id, reference_subset, disk_usage,
                                      disk=disk_usage, cores=8, memory='16G')
                filtered_bam_file_ids[dtype][reference_subset].append(j.rv())

    # IsoSeq hints
    iso_seq_hints_file_ids = []
    if len(iso_seq_file_ids) > 0:
        for bam_file_id, bai_file_id in iso_seq_file_ids:
            disk_usage = tools.toilInterface.find_total_disk_usage([bam_file_id, bai_file_id])
            j = job.addChildJobFn(generate_iso_seq_hints, bam_file_id, bai_file_id, disk=disk_usage)
            iso_seq_hints_file_ids.append(j.rv())

    # annotation hints
    if input_file_ids['annotation'] is not None:
        disk_usage = tools.toilInterface.find_total_disk_usage(input_file_ids['annotation'])
        j = job.addChildJobFn(generate_annotation_hints, input_file_ids['annotation'], disk=disk_usage)
        annotation_hints_file_id = j.rv()
    else:
        annotation_hints_file_id = None
    return job.addFollowOnJobFn(merge_bams, filtered_bam_file_ids, annotation_hints_file_id,
                                iso_seq_hints_file_ids).rv()


def namesort_bam(job, bam_file_id, bai_file_id, reference_subset, disk_usage, num_reads=50 ** 6):
    """
    Slices out the reference subset from a BAM, name sorts that subset, then chunks the resulting reads up for
    processing by filterBam.
    """

    def write_bam(r, ns_handle):
        """Write to the path, returns file ID"""
        outf = tools.fileOps.get_tmp_toil_file()
        outf_h = pysam.Samfile(outf, 'wb', template=ns_handle)
        for rec in r:
            outf_h.write(rec)
        outf_h.close()
        return job.fileStore.writeGlobalFile(outf)

    bam_path = job.fileStore.readGlobalFile(bam_file_id)
    is_paired = bam_is_paired(bam_path)
    job.fileStore.readGlobalFile(bai_file_id, bam_path + '.bai')
    name_sorted = tools.fileOps.get_tmp_toil_file(suffix='name_sorted.bam')
    cmd = [['samtools', 'view', '-b', bam_path] + list(reference_subset),
           ['sambamba', 'sort', '-t', '8', '-m', '15G', '-o', '/dev/stdout', '-n', '/dev/stdin']]
    tools.procOps.run_proc(cmd, stdout=name_sorted)
    ns_handle = pysam.Samfile(name_sorted)
    filtered_file_ids = []
    r = []
    for i, (qname, reads) in enumerate(itertools.groupby(ns_handle, lambda x: x.qname)):
        r.extend(list(reads))
        if i != 0 and i % num_reads == 0:
            file_id = write_bam(r, ns_handle)
            j = job.addChildJobFn(filter_bam, file_id, is_paired, disk='4G', memory='2G')
            filtered_file_ids.append(j.rv())
            r = []
    # do the last bin, if its non-empty
    if len(r) > 0:
        file_id = write_bam(r, ns_handle)
        j = job.addChildJobFn(filter_bam, file_id, is_paired, disk='4G', memory='2G')
        filtered_file_ids.append(j.rv())
    return job.addFollowOnJobFn(merge_filtered_bams, filtered_file_ids, disk=disk_usage, memory='16G').rv()


def filter_bam(job, file_id, is_paired):
    """
    Filters a name-sorted bam, returns a bam re-sorted by position
    """
    bam_path = job.fileStore.readGlobalFile(file_id)
    assert os.path.getsize(bam_path) > 0
    tmp_filtered = tools.fileOps.get_tmp_toil_file()
    filter_cmd = ['filterBam', '--uniq', '--in', bam_path, '--out', tmp_filtered]
    if is_paired is True:
        filter_cmd.extend(['--paired', '--pairwiseAlignments'])
    tools.procOps.run_proc(filter_cmd)
    assert os.path.getsize(tmp_filtered) > 0
    sort_tmp = tools.fileOps.get_tmp_toil_file()
    out_filter = tools.fileOps.get_tmp_toil_file()
    sort_cmd = ['samtools', 'sort', '-O', 'bam', '-T', sort_tmp, tmp_filtered]
    tools.procOps.run_proc(sort_cmd, stdout=out_filter)
    return job.fileStore.writeGlobalFile(out_filter)


def merge_filtered_bams(job, filtered_file_ids):
    """
    Merges filtered BAMs
    """
    local_paths = [job.fileStore.readGlobalFile(x) for x in filtered_file_ids]
    fofn = tools.fileOps.get_tmp_toil_file()
    with open(fofn, 'w') as outf:
        for l in local_paths:
            outf.write(l + '\n')
    out_bam = tools.fileOps.get_tmp_toil_file()
    cmd = ['samtools', 'merge', '-b', fofn, out_bam]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(out_bam)


def merge_bams(job, filtered_bam_file_ids, annotation_hints_file_id, iso_seq_hints_file_ids):
    """
    Takes a dictionary mapping reference chunks to filtered BAMs. For each reference chunk, these BAMs will be
    first concatenated then sorted, then passed off to hint building.
    """
    merged_bam_file_ids = {'BAM': {}, 'INTRONBAM': {}}
    for dtype in filtered_bam_file_ids:
        for ref_group, file_ids in filtered_bam_file_ids[dtype].iteritems():
            disk_usage = tools.toilInterface.find_total_disk_usage(file_ids)
            merged_bam_file_ids[dtype][ref_group] = job.addChildJobFn(cat_sort_bams, file_ids, disk=disk_usage).rv()

    return job.addFollowOnJobFn(build_hints, merged_bam_file_ids, annotation_hints_file_id, iso_seq_hints_file_ids).rv()


def cat_sort_bams(job, bam_file_ids):
    """
    Takes a list of bam file IDs and combines/sorts them.

    TODO: the 4096 file hack below is hacky. Should only be a problem for very fragmented references.
    """
    bamfiles = [job.fileStore.readGlobalFile(x) for x in bam_file_ids]
    # cat only 4095 bams at a time to avoid bash command length problems
    catfile = tools.fileOps.get_tmp_toil_file()
    sam_iter = tools.dataOps.grouper(bamfiles, 4095)
    # do the first one
    cmd = ['samtools', 'cat', '-o', catfile]
    cmd.extend(sam_iter.next())
    tools.procOps.run_proc(cmd)
    # do any subsequent ones left, creating a new file each time
    for more in sam_iter:
        old_catfile = catfile
        catfile = tools.fileOps.get_tmp_toil_file()
        cmd = ['samtools', 'cat', '-o', catfile, old_catfile]
        cmd.extend(more)
        tools.procOps.run_proc(cmd)
    # combine and merge
    merged = tools.fileOps.get_tmp_toil_file()
    sort_tmp = tools.fileOps.get_tmp_toil_file()
    cmd = ['samtools', 'sort', '-O', 'bam', '-T', sort_tmp, catfile]
    tools.procOps.run_proc(cmd, stdout=merged)
    return job.fileStore.writeGlobalFile(merged)


def build_hints(job, merged_bam_file_ids, annotation_hints_file_id, iso_seq_hints_file_ids):
    """
    Takes the merged BAM for a genome and produces both intron and exon hints.
    """
    intron_hints_file_ids = []
    exon_hints_file_ids = []
    for dtype in merged_bam_file_ids:
        for ref_group, file_ids in merged_bam_file_ids[dtype].iteritems():
            intron_hints_file_ids.append(job.addChildJobFn(build_intron_hints, file_ids).rv())
            if dtype == 'BAM':
                exon_hints_file_ids.append(job.addChildJobFn(build_exon_hints, file_ids).rv())
    disk_usage = tools.toilInterface.find_total_disk_usage(itertools.chain.from_iterable([intron_hints_file_ids,
                                                                                          exon_hints_file_ids,
                                                                                          [annotation_hints_file_id]]))
    return job.addFollowOnJobFn(cat_hints, intron_hints_file_ids, exon_hints_file_ids, annotation_hints_file_id,
                                iso_seq_hints_file_ids, disk=disk_usage).rv()


def build_intron_hints(job, merged_bam_file_id):
    """Builds intronhints from a BAM. Returns a fileID to the hints."""
    bam_file = job.fileStore.readGlobalFile(merged_bam_file_id)
    intron_gff_path = tools.fileOps.get_tmp_toil_file()
    cmd = ['bam2hints', '--intronsonly', '--in', bam_file, '--out', intron_gff_path]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(intron_gff_path)


def build_exon_hints(job, merged_bam_file_id):
    """Builds exonhints from a BAM Returns a fileID to the hints."""
    bam_file = job.fileStore.readGlobalFile(merged_bam_file_id)
    cmd = [['bam2wig', bam_file],
           ['wig2hints.pl', '--width=10', '--margin=10', '--minthresh=2', '--minscore=4', '--prune=0.1', '--src=W',
            '--type=ep', '--UCSC=/dev/null', '--radius=4.5', '--pri=4', '--strand=.']]
    exon_gff_path = tools.fileOps.get_tmp_toil_file()
    tools.procOps.run_proc(cmd, stdout=exon_gff_path)
    return job.fileStore.writeGlobalFile(exon_gff_path)


def generate_iso_seq_hints(job, bam_file_id, bai_file_id):
    """
    Generates hints from a IsoSeq BAM. Due to the usual depth of IsoSeq, there is no real need to split it up by
    chunks of reference sequence.

    Adapted from http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Augustus.PacBioGMAP
    """
    bam_path = job.fileStore.readGlobalFile(bam_file_id)
    job.fileStore.readGlobalFile(bai_file_id, bam_path + '.bai')
    pacbio_gff_path = tools.fileOps.get_tmp_toil_file()
    cmd = [['samtools', 'view', '-b', '-F', '4', bam_path],  # unmapped reads causes bamToPsl to crash
           ['bamToPsl', '/dev/stdin', '/dev/stdout'],
           ['sort', '-n', '-k', '16,16'],
           ['sort', '-s', '-k', '14,14'],
           ['perl', '-ne', '@f=split; print if ($f[0]>=100)'],
           ['blat2hints.pl', '--source=PB', '--nomult', '--ep_cutoff=20', '--in=/dev/stdin',
            '--out={}'.format(pacbio_gff_path)]]
    tools.procOps.run_proc(cmd)
    return job.fileStore.writeGlobalFile(pacbio_gff_path)


def generate_annotation_hints(job, annotation_hints_file_id):
    """
    Converts the annotation file into hints. First converts the gff3 directly to genePred so we can make use
    of the transcript library.

    Hints are derived from both CDS exonic intervals and intron intervals
    """
    annotation_gff3 = job.fileStore.readGlobalFile(annotation_hints_file_id)
    tm_gp = tools.fileOps.get_tmp_toil_file()
    cmd = ['gff3ToGenePred', '-rnaNameAttr=transcript_id', '-geneNameAttr=gene_id', '-honorStartStopCodons',
           annotation_gff3, tm_gp]
    tools.procOps.run_proc(cmd)
    tx_dict = tools.transcripts.get_gene_pred_dict(tm_gp)
    hints = []
    for tx_id, tx in tx_dict.iteritems():
        if tx.cds_size == 0:
            continue
        # rather than try to re-do the arithmetic, we will use the get_bed() function to convert this transcript
        cds_tx = tools.transcripts.Transcript(tx.get_bed(new_start=tx.thick_start, new_stop=tx.thick_stop))
        for intron in cds_tx.intron_intervals:
            r = [intron.chromosome, 'a2h', 'intron', intron.start + 1, intron.stop, 0, intron.strand, '.',
                 'grp={};src=M;pri=2'.format(tx_id)]
            hints.append(r)
        for exon in cds_tx.exon_intervals:
            r = [exon.chromosome, 'a2h', 'CDS', exon.start + 1, exon.stop, 0, exon.strand, '.',
                 'grp={};src=M;pri=2'.format(tx_id)]
            hints.append(r)
    annotation_hints_gff = tools.fileOps.get_tmp_toil_file()
    tools.fileOps.print_rows(annotation_hints_gff, hints)
    return job.fileStore.writeGlobalFile(annotation_hints_gff)


def cat_hints(job, intron_hints_file_ids, exon_hints_file_ids, annotation_hints_file_id, iso_seq_hints_file_ids):
    """Returns file ID to combined, sorted hints"""
    cat_hints = tools.fileOps.get_tmp_toil_file()
    with open(cat_hints, 'w') as outf:
        for file_id in itertools.chain(intron_hints_file_ids, exon_hints_file_ids):
            f = job.fileStore.readGlobalFile(file_id)
            for line in open(f):
                outf.write(line)
        if annotation_hints_file_id is not None:
            f = job.fileStore.readGlobalFile(annotation_hints_file_id)
            for line in open(f):
                outf.write(line)
    # sorted so that hints that should be summarized are below each other
    cmd = [['sort', '-n', '-k4,4', cat_hints],
           ['sort', '-s', '-n', '-k5,5'],
           ['sort', '-s', '-k3,3'],
           ['sort', '-s', '-k1,1'],
           ['join_mult_hints.pl']]
    combined_hints = tools.fileOps.get_tmp_toil_file()
    tools.procOps.run_proc(cmd, stdout=combined_hints)
    # don't add the IsoSeq until after join_mult_hints because we don't want them to be joined
    with open(combined_hints, 'a') as outf:
        for file_id in iso_seq_hints_file_ids:
            f = job.fileStore.readGlobalFile(file_id)
            for line in open(f):
                outf.write(line)
    # sort the combined hints, now sorting by chrom and start
    sorted_combined_hints = tools.fileOps.get_tmp_toil_file()
    tools.misc.sort_gff(combined_hints, sorted_combined_hints)
    return job.fileStore.writeGlobalFile(sorted_combined_hints)


###
# Functions
###


def validate_bam_fasta_pairs(bam_path, fasta_sequences, genome):
    """
    Make sure that this BAM is actually aligned to this fasta. Every sequence should be the same length. Sequences
    can exist in the reference that do not exist in the BAM, but not the other way around.
    """
    handle = pysam.Samfile(bam_path, 'rb')
    bam_sequences = {(n, s) for n, s in zip(*[handle.references, handle.lengths])}
    difference = bam_sequences - fasta_sequences
    if len(difference) > 0:
        base_err = 'Error: BAM {} has the following sequence/length pairs not found in the {} fasta: {}.'
        err = base_err.format(bam_path, genome, ','.join(['-'.join(map(str, x)) for x in difference]))
        raise UserException(err)
    missing_seqs = fasta_sequences - bam_sequences
    if len(missing_seqs) > 0:
        base_msg = 'BAM {} does not have the following sequence/length pairs in its header: {}.'
        msg = base_msg.format(bam_path, ','.join(['-'.join(map(str, x)) for x in missing_seqs]))
        logger.warning(msg)


def bam_is_paired(bam_path, num_reads=20000, paired_cutoff=0.75):
    """
    Infers the paired-ness of a bam file.
    """
    sam = pysam.Samfile(bam_path)
    count = 0
    for rec in itertools.islice(sam, num_reads):
        if rec.is_paired:
            count += 1
    if tools.mathOps.format_ratio(count, num_reads) > 0.75:
        return True
    elif tools.mathOps.format_ratio(count, num_reads) < 1 - paired_cutoff:
        return False
    else:
        raise UserException("Unable to infer pairing from bamfile {}".format(bam_path))


def group_references(sam_handle, num_bases=10 ** 7, max_seqs=1000):
    """
    Group up references by num_bases, unless that exceeds max_seqs. A greedy implementation of the bin packing problem.
    """
    name_iter = itertools.izip(*[sam_handle.references, sam_handle.lengths])
    name, size = name_iter.next()
    this_bin = [name]
    bin_base_count = size
    num_seqs = 1
    for name, size in name_iter:
        bin_base_count += size
        num_seqs += 1
        if bin_base_count >= num_bases or num_seqs > max_seqs:
            yield this_bin
            this_bin = [name]
            bin_base_count = size
            num_seqs = 1
        else:
            this_bin.append(name)
    yield this_bin
