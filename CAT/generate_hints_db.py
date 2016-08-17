"""
Generate a hints database file from RNAseq alignments for AugustusTMR/AugustusCGP.

Expects a config file to be passed in with the paths to the files. Example:

[FASTA]
genome1 = /path/to/fasta

[ANNOTATION]
annotation = /path/to/gff3

[BAM]
genome1 = /path/to/bam1.bam, /path/to/bam2.bam

The annotation field is optional, but will help AugustusCGP make better predictions.

"""
import os
import shutil
import luigi
import pysam
import itertools
import tempfile
from configobj import ConfigObj
from toil.job import Job
from toil.common import Toil

import tools.toilInterface
import tools.fileOps
import tools.procOps
import tools.mathOps
import tools.transcripts
import luigi.contrib.sqla
from luigi.util import inherits, requires


class UserException(Exception):
    pass


class BuildHints(luigi.WrapperTask):
    """
    Main entry point. Parses input files, and launches the next task.
    """
    # path to config file
    config = luigi.Parameter(default='hints_cfg.cfg')
    augustus_hints_db = luigi.Parameter(default='augustus_hints.db')
    work_dir = luigi.Parameter(default='hints_work')
    # toil parameters, which match tools.toilInterface.ToilTask
    workDir = luigi.Parameter(default=tempfile.gettempdir(), significant=False)
    batchSystem = luigi.Parameter(default='singleMachine', significant=False)
    maxCores = luigi.IntParameter(default=16, significant=False)
    logLevel = luigi.Parameter(default='WARNING', significant=False)
    cleanWorkDir = luigi.Parameter(default=None, significant=False)  # debugging option

    def parse_cfg(self):
        # configspec validates the input config file
        configspec = ['[FASTA]', '__many__ = string', '[ANNOTATION]', '__many__ = string', '[BAM]', '__many__ = list']
        parser = ConfigObj(self.config, configspec=configspec)
        # if a given genome only has one BAM, it is a string. Fix this.
        for genome in parser['BAM']:
            if isinstance(parser['BAM'][genome], str):
                parser['BAM'][genome] = [parser['BAM'][genome]]
        # do some input validation
        for genome in parser['BAM']:
            for bam in parser['BAM'][genome]:
                assert os.path.exists(bam)
                assert os.path.exists(bam + '.bai')
        # TODO: validate all inputs. Some validation happens in the toil pipeline as well.
        return parser

    def requires(self):
        cfg = self.parse_cfg()
        hint_paths = {}
        flat_fasta_paths = {}
        for genome in cfg['FASTA']:
            flat_fasta = self.clone(GenomeFlatFasta, genome=genome, cfg=cfg)
            yield flat_fasta
            flat_fasta_paths[genome] = flat_fasta.output().path
            annotation = cfg['ANNOTATION'].get(genome, None)
            hints = self.clone(GenerateHints, genome=genome, flat_fasta=flat_fasta.output().path,
                               annotation=annotation, cfg=cfg)
            yield hints
            hint_paths[genome] = hints.output().path
        yield self.clone(BuildDb, cfg=cfg, hint_paths=hint_paths, flat_fasta_paths=flat_fasta_paths)


@inherits(BuildHints)
class GenomeFlatFasta(luigi.Task):
    """
    Flattens a genome fasta using pyfasta, copying it to the work directory. Requires the pyfasta package.
    """
    genome = luigi.Parameter()
    cfg = luigi.Parameter()

    def output(self):
        path = os.path.abspath(os.path.join(self.work_dir, self.genome + '.fasta'))
        tools.fileOps.ensure_file_dir(path)
        return luigi.LocalTarget(path)

    def run(self):
        shutil.copyfile(self.cfg['FASTA'][self.genome], self.output().path)
        cmd = ['pyfasta', 'flatten', self.output().path]
        tools.procOps.run_proc(cmd)


@requires(GenomeFlatFasta)
class GenerateHints(tools.toilInterface.ToilTask):
    """
    Generate hints for each genome as a separate Toil pipeline.
    """
    genome = luigi.Parameter()
    flat_fasta = luigi.Parameter()
    annotation = luigi.Parameter()
    cfg = luigi.Parameter()

    def output(self):
        hints = os.path.abspath(os.path.join(self.work_dir, self.genome + '.reduced_hints.gff'))
        return luigi.LocalTarget(hints)

    def run(self):
        job_store = os.path.join(self.work_dir, 'toil', self.genome)
        toil_options = self.prepare_toil_options(job_store)
        toil_options.defaultMemory = '8G'  # TODO: don't hardcode this.
        completed = generate_hints(self.genome, self.flat_fasta, self.cfg, self.annotation, self.output().path,
                                   toil_options)
        if completed is False:  # we did not have any hints to generate for this genome
            self.output().open('w').close()


@inherits(BuildHints)
class BuildDb(luigi.Task):
    """
    Constructs the hints database from a series of reduced hints GFFs

    TODO: output() should be way smarter than this. Stefanie is going to work on flags that modify the load2sqlite
    program
    """
    cfg = luigi.Parameter()
    hint_paths = luigi.Parameter()
    flat_fasta_paths = luigi.Parameter()

    def requires(self):
        for genome in self.cfg['FASTA']:
            flat_fasta = self.clone(GenomeFlatFasta, genome=genome, cfg=self.cfg)
            annotation = self.cfg['ANNOTATION'].get(genome, None)
            hints = self.clone(GenerateHints, genome=genome, flat_fasta=flat_fasta.output().path,
                               annotation=annotation, cfg=self.cfg)
            yield hints

    def output(self):
        tools.fileOps.ensure_file_dir(self.augustus_hints_db)
        # check first if each genome is in speciesnames
        # then check if # of hints for the genome matches the # of hints in the gff
        return luigi.LocalTarget('/thisfiledoesnotexist')  # TODO: this will just always run now.

    def run(self):
        for genome in self.cfg['FASTA']:
            base_cmd = ['load2sqlitedb', '--noIdx', '--clean', '--species={}'.format(genome),
                        '--dbaccess={}'.format(self.augustus_hints_db)]
            tools.procOps.run_proc(base_cmd + [self.flat_fasta_paths[genome]])
            if genome in self.hint_paths and os.path.getsize(self.hint_paths[genome]) > 0:
                tools.procOps.run_proc(base_cmd + [self.hint_paths[genome]])
        cmd = ['load2sqlitedb', '--makeIdx', '--clean', '--dbaccess={}'.format(self.augustus_hints_db)]
        tools.procOps.run_proc(cmd)


###
# Toil pipeline
###


def generate_hints(genome, flat_fasta, cfg, annotation, out_gff_path, toil_options):
    """
    Entry point for hints database Toil pipeline.
    """
    with Toil(toil_options) as toil:
        if not toil.options.restart:
            fasta_file_id = toil.importFile('file://' + flat_fasta)
            fasta_gdx_file_id = toil.importFile('file://' + flat_fasta + '.gdx')
            fasta_flat_file_id = toil.importFile('file://' + flat_fasta + '.flat')
            if genome in cfg['BAM'] and cfg['BAM'][genome] is not None:
                bam_file_ids = {}
                for bam_path in cfg['BAM'][genome]:
                    bam_file_ids[os.path.basename(bam_path)] = (toil.importFile('file://' + bam_path),
                                                                toil.importFile('file://' + bam_path + '.bai'))
            else:
                bam_file_ids = None
            annotation_file_id = toil.importFile('file://' + annotation) if annotation is not None else None
            input_file_ids = {'genome_fasta': fasta_file_id, 'genome_gdx': fasta_gdx_file_id,
                              'genome_flat': fasta_flat_file_id,
                              'bams': bam_file_ids, 'annotation': annotation_file_id}
            job = Job.wrapJobFn(setup_hints, genome, input_file_ids)
            combined_hints = toil.start(job)
        else:
            combined_hints = toil.restart()
        if combined_hints is not None:
            tools.fileOps.ensure_file_dir(out_gff_path)
            toil.exportFile(combined_hints, 'file://' + out_gff_path)
            return True
        else:
            return False


def setup_hints(job, genome, input_file_ids):
    """
    Generates hints for a given genome with a list of BAMs. Will add annotation if it exists.
    """
    # import FASTA, flatten it, export flattened fasta
    genome_fasta = tools.toilInterface.load_fasta_from_filestore(job, input_file_ids['genome_fasta'],
                                                                 input_file_ids['genome_gdx'],
                                                                 input_file_ids['genome_flat'],
                                                                 prefix='genome', upper=False)
    if input_file_ids['bams'] is not None:
        # validate BAMs
        fasta_sequences = {x.split()[0] for x in genome_fasta.keys()}
        filtered_bam_file_ids = []
        for original_path, (bam_file_id, bai_file_id) in input_file_ids['bams'].iteritems():
            bam_path = job.fileStore.readGlobalFile(bam_file_id)
            validate_bam_fasta_pairs(original_path, bam_path, fasta_sequences, genome)
            sam_handle = pysam.Samfile(bam_path)
            paired = ['--paired', '--pairwiseAlignments'] if bam_is_paired(original_path, bam_path) is True else []
            for reference_subset in group_references(sam_handle):
                j = job.addChildJobFn(sort_by_name, bam_file_id, bai_file_id, reference_subset, paired)
                filtered_bam_file_ids.append(j.rv())
    else:
        filtered_bam_file_ids = []
    if input_file_ids['annotation'] is not None:
        j = job.addChildJobFn(generate_annotation_hints, input_file_ids['annotation'])
        annotation_hints = j.rv()
    else:
        annotation_hints = None
    # returns path to filtered GFF output
    return job.addFollowOnJobFn(merge_bams, filtered_bam_file_ids, annotation_hints).rv()


def sort_by_name(job, bam_file_id, bai_file_id, reference_subset, paired):
    """
    Slices out a chromosome from a BAM, re-sorts by name, filters the reads, then re-sorts by position.
    filterBam does weird things when piped to stdout, so I don't do that.
    """
    tmp_filtered = tools.fileOps.get_tmp_toil_file(suffix='filtered.bam')
    bam_path = job.fileStore.readGlobalFile(bam_file_id)
    job.fileStore.readGlobalFile(bai_file_id, bam_path + '.bai')
    sort_tmp = tools.fileOps.get_tmp_toil_file()
    cmd = [['samtools', 'view', '-b', bam_path],
           ['samtools', 'sort', '-O', 'bam', '-T', sort_tmp, '-n', '-l', '0', '-'],
           ['filterBam', '--uniq', '--in', '/dev/stdin', '--out', tmp_filtered]]
    cmd[-1].extend(paired)
    cmd[0].extend(reference_subset)
    tools.procOps.run_proc(cmd)
    out_filter = tools.fileOps.get_tmp_toil_file(suffix='sorted.filtered.bam')
    sort_cmd = ['samtools', 'sort', '-O', 'bam', '-T', sort_tmp, tmp_filtered]
    tools.procOps.run_proc(sort_cmd, stdout=out_filter)
    tools.procOps.run_proc(['samtools', 'index', out_filter])
    filtered_bam_file_id = job.fileStore.writeGlobalFile(out_filter)
    filtered_bam_index_file_id = job.fileStore.writeGlobalFile(out_filter + '.bai')
    return filtered_bam_file_id, filtered_bam_index_file_id


def merge_bams(job, filtered_bam_file_ids, annotation_hints):
    """
    Takes the filtered BAMs that were split up by sequences and merges them into one mega-BAM.
    TODO: this could probably be made faster through a hierarchical merging process. Not sure how samtools merge
    sorts in the back end.
    TODO: this should probably request a sizeable amount of HDD space, that scales based on input data.
    """
    if len(filtered_bam_file_ids) > 0:
        # import every single BAM
        bam_paths = []
        for filtered_bam_file_id, filtered_bam_index_file_id in filtered_bam_file_ids:
            bam = tools.fileOps.get_tmp_toil_file()
            bai = bam + '.bai'
            job.fileStore.readGlobalFile(filtered_bam_file_id, bam)
            job.fileStore.readGlobalFile(filtered_bam_index_file_id, bai)
            bam_paths.append(bam)
        # write a fofn
        fofn = tools.fileOps.get_tmp_toil_file()
        tools.fileOps.print_rows(fofn, bam_paths, sep='')
        # generate output
        merged_path = tools.fileOps.get_tmp_toil_file(suffix='bam')
        cmd = ['samtools', 'merge', '-b', fofn, merged_path]
        tools.procOps.run_proc(cmd)
        merged_bam_file_id = job.fileStore.writeGlobalFile(merged_path)
    else:
        merged_bam_file_id = None
    return job.addFollowOnJobFn(build_hints, merged_bam_file_id, annotation_hints).rv()


def build_hints(job, merged_bam_file_id, annotation_hints):
    """
    Takes the merged BAM for a genome and produces both intron and exon hints.
    """
    intron_hints = job.addChildJobFn(build_intron_hints, merged_bam_file_id).rv()
    exon_hints = job.addChildJobFn(build_exon_hints, merged_bam_file_id).rv()
    return job.addFollowOnJobFn(cat_hints, intron_hints, exon_hints, annotation_hints).rv()


def build_intron_hints(job, merged_bam_file_id):
    """Builds intronhints from a BAM. Returns a fileID to the hints."""
    if merged_bam_file_id is not None:
        bam_file = job.fileStore.readGlobalFile(merged_bam_file_id)
        intron_gff_path = tools.fileOps.get_tmp_toil_file(suffix='intron_hints.gff')
        cmd = ['bam2hints', '--intronsonly', '--in', bam_file, '--out', intron_gff_path]
        tools.procOps.run_proc(cmd)
        return job.fileStore.writeGlobalFile(intron_gff_path)
    else:
        return None


def build_exon_hints(job, merged_bam_file_id):
    """Builds exonhints from a BAM Returns a fileID to the hints."""
    if merged_bam_file_id is not None:
        bam_file = job.fileStore.readGlobalFile(merged_bam_file_id)
        cmd = [['bam2wig', bam_file],
               ['wig2hints.pl', '--width=10', '--margin=10', '--minthresh=2', '--minscore=4', '--prune=0.1', '--src=W',
                '--type=ep', '--UCSC=/dev/null', '--radius=4.5', '--pri=4', '--strand="."']]
        exon_gff_path = tools.fileOps.get_tmp_toil_file(suffix='exon_hints.gff')
        tools.procOps.run_proc(cmd, stdout=exon_gff_path)
        return job.fileStore.writeGlobalFile(exon_gff_path)
    else:
        return None


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


def cat_hints(job, intron_hints, exon_hints, annotation_hints):
    """Returns file ID to combined, sorted hints"""
    if all(x is None for x in (intron_hints, exon_hints, annotation_hints)):
        return None  # we are just loading the genome
    cmd = [['cat'],
           ['sort', '-n', '-k4,4'],
           ['sort', '-s', '-n', '-k5,5'],
           ['sort', '-s', '-n', '-k3,3'],
           ['sort', '-s', '-k1,1'],
           ['join_mult_hints.pl']]
    for hint_file_id in [intron_hints, exon_hints, annotation_hints]:
        if hint_file_id is not None:
            hints_path = job.fileStore.readGlobalFile(hint_file_id)
            cmd[0].append(hints_path)
    combined_hints = tools.fileOps.get_tmp_toil_file()
    tools.procOps.run_proc(cmd, stdout=combined_hints)
    return job.fileStore.writeGlobalFile(combined_hints)


###
# Functions
###


def validate_bam_fasta_pairs(original_path, bam_path, fasta_sequences, genome):
    """
    Make sure that this BAM is actually aligned to this fasta
    """
    handle = pysam.Samfile(bam_path, 'rb')
    if set(handle.references) != fasta_sequences:
        base_err = 'Error: BAM {} does not have the same sequences as the FASTA for genome {}'
        err = base_err.format(original_path, genome)
        raise UserException(err)


def bam_is_paired(original_path, bam_path, num_reads=100000, paired_cutoff=0.75):
    """
    Infers the paired-ness of a bam file.
    """
    sam = pysam.Samfile(bam_path)
    count = 0
    for i, rec in enumerate(sam):
        if rec.is_paired:
            count += 1
        if i == num_reads:
            break
    if tools.mathOps.format_ratio(count, num_reads) > 0.75:
        return True
    elif tools.mathOps.format_ratio(count, num_reads) < 1 - paired_cutoff:
        return False
    else:
        raise UserException("Unable to infer pairing from bamfile {}".format(original_path))


def group_references(sam_handle, num_bases=10 ** 6, max_seqs=100):
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
