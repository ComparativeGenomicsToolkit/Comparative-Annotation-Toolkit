#!/usr/bin/env python2.7
"""
Comparative Annotation Toolkit.
"""
import itertools
import os
import luigi
import tempfile
import pandas as pd
import tools.fileOps
import tools.procOps
import tools.toilInterface
import tools.hal
import tools.transcripts
import tools.psl
import tools.gff3
from tools.luigiAddons import multiple_inherits, multiple_requires
from luigi.util import inherits, requires
from src.abstractClasses import AbstractAtomicFileTask
from src.chaining import chaining


class UserException(Exception):
    pass


class RunCat(luigi.WrapperTask):
    """
    Task that executes the entire pipeline.
    """
    hal = luigi.Parameter()
    ref_genome = luigi.Parameter()
    annotation = luigi.Parameter()
    out_dir = luigi.Parameter(default='./cat_output')
    work_dir = luigi.Parameter(default=tempfile.gettempdir())
    target_genomes = luigi.TupleParameter(default=None)
    augustus = luigi.BoolParameter(default=False)
    augustus_genomes = luigi.TupleParameter(default=None)
    augustus_hints_db = luigi.Parameter(default=None)

    @staticmethod
    def resolve_target_genomes(hal, target_genomes):
        if target_genomes is None:
            target_genomes = tools.hal.extract_genomes(hal)
        return target_genomes

    def requires(self):
        yield self.clone(PrepareFiles)
        yield self.clone(Chaining)
        yield self.clone(TransMap)
        #yield EvalTransMap()
        #yield AugustusTMR()
        #yield Align()
        #yield EvalTranscripts()
        #yield Consensus()
        #yield Plots()


@inherits(RunCat)
class PrepareFiles(luigi.WrapperTask):
    """
    Wrapper for file preparation tasks GenomeFiles and ReferenceFiles
    """
    def validate(self):
        # TODO: make sure that all input args exist
        pass

    def requires(self):
        self.validate()
        yield self.clone(GenomeFiles)
        yield self.clone(ReferenceFiles)


@inherits(PrepareFiles)
class GenomeFiles(luigi.WrapperTask):
    """
    WrapperTask for producing all genome files.

    GenomeFiles -> GenomeFasta -> GenomeTwoBit -> GenomeFlatFasta -> GenomeFastaIndex
                -> GenomeSizes

    """
    @staticmethod
    def get_args(work_dir, genome):
        base_dir = os.path.join(work_dir, 'genome_files')
        args = {'genome': genome,
                 'fasta': os.path.join(base_dir, genome + '.fa'),
                 'two_bit': os.path.join(base_dir, genome + '.2bit'),
                 'fasta_index': os.path.join(base_dir, genome + '.fa.fai'),
                 'sizes': os.path.join(base_dir, genome + '.chrom.sizes'),
                 'flat_fasta': os.path.join(base_dir, genome + '.fa.flat')}
        return args

    def requires(self):
        target_genomes = RunCat.resolve_target_genomes(self.hal, self.target_genomes)
        for genome in itertools.chain(target_genomes, [self.ref_genome]):
            args = self.get_args(self.work_dir, genome)
            yield self.clone(GenomeFastaIndex, **args)
            yield self.clone(GenomeSizes, **args)


@inherits(PrepareFiles)
class GenomeFasta(AbstractAtomicFileTask):
    """
    Produce a fasta file from a hal file. Requires hal2fasta.
    """
    genome = luigi.Parameter()
    fasta = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.fasta)

    def run(self):
        cmd = ['hal2fasta', self.hal, self.genome]
        self.run_cmd(cmd)


@requires(GenomeFasta)
class GenomeTwoBit(AbstractAtomicFileTask):
    """
    Produce a 2bit file from a fasta file. Requires kent tool faToTwoBit.
    Needs to be done BEFORE we flatten.
    """
    two_bit = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.two_bit)

    def run(self):
        cmd = ['faToTwoBit', self.fasta, '/dev/stdout']
        self.run_cmd(cmd)


@inherits(PrepareFiles)
class GenomeSizes(AbstractAtomicFileTask):
    """
    Produces a genome chromosome sizes file. Requires halStats.
    """
    genome = luigi.Parameter()
    sizes = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.sizes)

    def run(self):
        cmd = ['halStats', '--chromSizes', self.genome, self.hal]
        self.run_cmd(cmd)


@requires(GenomeTwoBit)
class GenomeFlatFasta(AbstractAtomicFileTask):
    """
    Flattens a genome fasta in-place using pyfasta. Requires the pyfasta package.
    """
    flat_fasta = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.flat_fasta)

    def run(self):
        cmd = ['pyfasta', 'flatten', self.fasta]
        tools.procOps.run_proc(cmd)


@requires(GenomeFlatFasta)
class GenomeFastaIndex(AbstractAtomicFileTask):
    """
    Index a flat fasta. This must be done after flattening.
    """
    fasta_index = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.fasta_index)

    def run(self):
        cmd = ['samtools', 'faidx', self.fasta]
        tools.procOps.run_proc(cmd)


@inherits(PrepareFiles)
class ReferenceFiles(luigi.WrapperTask):
    """
    WrapperTask for producing annotation files.

    ReferenceFiles -> Gff3ToGenePred -> TranscriptBed -> TranscriptFasta -> FlatTranscriptFasta
                            V
                         FakePsl

    """
    @staticmethod
    def get_args(work_dir, annotation):
        base_dir = os.path.join(work_dir, 'reference')
        args = {'annotation_gene_pred': os.path.join(base_dir, annotation + '.gp'),
                 'annotation_attrs': os.path.join(base_dir, annotation + '.attrs.tsv'),
                 'transcript_fasta':  os.path.join(base_dir, annotation + '.fa'),
                 'transcript_flat_fasta': os.path.join(base_dir, annotation + '.fa.flat'),
                 'transcript_bed': os.path.join(base_dir, annotation + '.bed'),
                 'fake_psl': os.path.join(base_dir, annotation + '.psl')}
        return args

    def requires(self):
        args = self.get_args(self.work_dir, self.annotation)
        args.update(GenomeFiles.get_args(self.work_dir, self.ref_genome))
        yield self.clone(FakePsl, **args)
        yield self.clone(FlatTranscriptFasta, **args)
        yield self.clone(Gff3ToAttrs, **args)


@inherits(ReferenceFiles)
class Gff3ToGenePred(luigi.Task):
    """
    Generates a genePred from a gff3 file. Also produces an attributes table.
    TODO: I really should just create the genePred from the gff3 parser yourself, instead of this hack.
    However, then I will spend a week writing my own converter and probably get it wrong.
    """
    annotation_gene_pred = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.annotation_gene_pred)

    def run(self):
        tmp_gp = luigi.LocalTarget(is_tmp=True)
        cmd = ['gff3ToGenePred', '-honorStartStopCodons', self.annotation, tmp_gp.path]
        tools.procOps.run_proc(cmd)
        fixed_gp = self.munge_gp(tmp_gp)
        tools.fileOps.atomic_install(fixed_gp.path, self.annotation_gene_pred)


@inherits(ReferenceFiles)
class Gff3ToAttrs(luigi.Task):
    """
    Uses the gff3 parser to extract the attributes table.
    """
    annotation_gene_pred = luigi.Parameter()
    annotation_attrs = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.annotation_attrs)

    def run(self):
        results = tools.gff3.extract_attrs(self.annotation)
        with self.output().open('w') as outf:
            tools.fileOps.print_rows(outf, results.itervalues())


@requires(Gff3ToGenePred)
class TranscriptBed(AbstractAtomicFileTask):
    """
    Produces a BED record from the input genePred annotation. Makes use of Kent tool genePredToBed
    """
    transcript_bed = luigi.Parameter()
    annotation_gene_pred = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.transcript_bed)

    def run(self):
        cmd = ['genePredToBed', self.annotation_gene_pred, '/dev/stdout']
        self.run_cmd(cmd)


@multiple_requires(GenomeFlatFasta, TranscriptBed)
class TranscriptFasta(AbstractAtomicFileTask):
    """
    Produces a fasta for each transcript. Requires bedtools.
    """
    transcript_fasta = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.transcript_fasta)

    def run(self):
        cmd = ['bedtools', 'getfasta', '-fi', self.fasta, '-bed', self.transcript_bed, '-fo', '/dev/stdout',
               '-name', '-split', '-s']
        self.run_cmd(cmd)


@requires(TranscriptFasta)
class FlatTranscriptFasta(AbstractAtomicFileTask):
    """
    Flattens the transcript fasta for pyfasta.
    """
    transcript_fasta = luigi.Parameter()
    transcript_flat_fasta = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.transcript_flat_fasta)

    def run(self):
        cmd = ['pyfasta', 'flatten', self.transcript_fasta]
        tools.procOps.run_proc(cmd)


@multiple_requires(Gff3ToGenePred, GenomeSizes)
class FakePsl(AbstractAtomicFileTask):
    """
    Produces a fake PSL mapping transcripts to the genome, using the Kent tool genePredToFakePsl
    """
    fake_psl = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.fake_psl)

    def run(self):
        cmd = ['genePredToFakePsl', '-chromSize={}'.format(self.sizes), 'noDB',
               self.annotation_gene_pred, '/dev/stdout', '/dev/null']
        self.run_cmd(cmd)


@inherits(RunCat)
class Chaining(luigi.WrapperTask):
    """
    WrapperTask for producing chain files for each combination of ref_genome-target_genome
    """
    @staticmethod
    def get_args(work_dir, ref_genome, genome):
        base_dir = os.path.join(work_dir, 'chaining')
        args = {'ref_genome': ref_genome, 'genome': genome,
                'chain_file': os.path.join(base_dir, '{}-{}.chain'.format(ref_genome, genome))}
        ref_files = GenomeFiles.get_args(work_dir, ref_genome)
        tgt_files = GenomeFiles.get_args(work_dir, genome)
        args['query_two_bit'] = ref_files['two_bit']
        args['target_two_bit'] = tgt_files['two_bit']
        args['query_sizes'] = ref_files['sizes']
        return args

    def validate(self):
        # TODO: make sure that all input args exist
        pass

    def requires(self):
        self.validate()
        target_genomes = RunCat.resolve_target_genomes(self.hal, self.target_genomes)
        for target_genome in target_genomes:
            args = self.get_args(self.work_dir, self.ref_genome, target_genome)
            yield self.clone(PairwiseChaining, chain_args=args, genome=target_genome)


@inherits(Chaining)
class PairwiseChaining(tools.toilInterface.ToilTask):
    """
    Executes Kent-style genomic chaining on a pair of genomes, producing a chain file.
    """
    chain_args = luigi.DictParameter()  # dict to pass directly to chaining toil module
    genome = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.chain_args['chain_file'])

    def requires(self):
        return self.clone(GenomeFiles)

    def run(self):
        job_sub_path = os.path.join('chaining', self.genome)
        toil_options = self.prepare_toil_options(job_sub_path)
        chaining(self.chain_args, self.hal, toil_options)


@inherits(RunCat)
class TransMap(luigi.WrapperTask):
    """
    Runs transMap.
    """
    @staticmethod
    def get_args(work_dir, genome, ref_genome, annotation):
        base_dir = os.path.join(work_dir, 'transMap')
        ref_files = ReferenceFiles.get_args(work_dir, annotation)
        tgt_genome_files = GenomeFiles.get_args(work_dir, genome)
        chain_args = Chaining.get_args(work_dir, ref_genome, genome)
        args = {'two_bit': tgt_genome_files['two_bit'],
                'chain_file': chain_args['chain_file'],
                'transcript_fasta': ref_files['transcript_fasta'], 'fake_psl': ref_files['fake_psl'],
                'annotation_gene_pred': ref_files['annotation_gene_pred'],
                'tm_psl': os.path.join(base_dir, genome + '.psl'), 'tm_gp': os.path.join(base_dir, genome + '.gp')}
        return args

    def validate(self):
        # TODO: make sure that all input args exist
        pass

    def requires(self):
        self.validate()
        target_genomes = RunCat.resolve_target_genomes(self.hal, self.target_genomes)
        for target_genome in target_genomes:
            args = self.get_args(self.work_dir, target_genome, self.ref_genome, self.annotation)
            yield self.clone(TransMapGp, tm_args=args, genome=target_genome)


@inherits(RunCat)
class TransMapPsl(luigi.Task):
    """
    Runs transMap. Requires Kent tools pslMap, postTransMapChain, pslRecalcMatch
    """
    tm_args = luigi.DictParameter()
    genome = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.tm_args['tm_psl'])

    def requires(self):
        return self.clone(PrepareFiles), self.clone(Chaining)

    def run(self):
        tools.fileOps.ensure_file_dir(self.output().path)
        psl_cmd = ['pslMap', '-chainMapFile', self.tm_args['fake_psl'],
                   self.tm_args['chain_file'], '/dev/stdout']
        post_chain_cmd = ['postTransMapChain', '/dev/stdin', '/dev/stdout']
        sort_cmd = ['sort', '-k', '14,14', '-k', '16,16n']
        recalc_cmd = ['pslRecalcMatch', '/dev/stdin', self.tm_args['two_bit'], self.tm_args['transcript_fasta'],
                      'stdout']
        cmd_list = [psl_cmd, post_chain_cmd, sort_cmd, recalc_cmd]
        # hacky way to make unique - capture output to a file, then process
        tmp_file = luigi.LocalTarget(is_tmp=True)
        with tmp_file.open('w') as tmp_fh:
            tools.procOps.run_proc(cmd_list, stdout=tmp_fh)
        with self.output().open('w') as outf:
            for q_name, psl_rec in tools.psl.psl_iterator(tmp_file.path, make_unique=True):
                outf.write('\t'.join(map(str, psl_rec.psl_string())) + '\n')


@requires(TransMapPsl)
class TransMapGp(AbstractAtomicFileTask):
    """
    Produces the final transMapped genePred
    """
    def output(self):
        return luigi.LocalTarget(self.tm_args['tm_gp'])

    def run(self):
        cmd = ['transMapPslToGenePred', '-nonCodingGapFillMax=50', '-codingGapFillMax=50',
               self.tm_args['annotation_gene_pred'], self.tm_args['tm_psl'], '/dev/stdout']
        self.run_cmd(cmd)


@inherits(RunCat)
class EvalTransMap(luigi.WrapperTask):
    """
    WrapperTask for first evaluation step. This step evaluates transMap alignments for having problems such as
    invalid splice sites, losing original introns, original start/stop, scaffold gaps
    """
    def validate(self):
        # TODO: make sure that all input args exist
        pass


if __name__ == '__main__':
    luigi.build([Chaining(hal='1509.hal',
                target_genomes=('C57BL_6NJ',), work_dir='test', ref_genome='C57B6J',
                annotation='Mus_musculus.GRCm38.83.gff3')], logging_conf_file='logging.cfg')

