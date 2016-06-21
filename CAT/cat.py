#!/usr/bin/env python2.7
"""
Comparative Annotation Toolkit.
"""
import itertools
import os
import luigi
import tempfile
import tools.fileOps
import tools.procOps
import tools.toilInterface
import tools.hal
import tools.transcripts
import tools.psl
import tools.gff3
import tools.misc
from tools.luigiAddons import multiple_inherits, multiple_requires, AbstractAtomicFileTask
from luigi.util import inherits, requires
from chaining import chaining
from augustus import augustus


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
    tm_cfg = luigi.Parameter(default='augustus_cfgs/extrinsic.ETM1.cfg', significant=False)
    tmr_cfg = luigi.Parameter(default='augustus_cfgs/extrinsic.ETM2.cfg', significant=False)

    @staticmethod
    def resolve_target_genomes(hal, target_genomes):
        if target_genomes is None:
            target_genomes = tools.hal.extract_genomes(hal)
        return target_genomes

    @staticmethod
    def resolve_augustus_genomes(hal, target_genomes, augustus_genomes):
        """Augustus genomes can either be A) all genomes, B) target genomes, C) specified separately"""
        if augustus_genomes is not None:
            return augustus_genomes
        return RunCat.resolve_target_genomes(hal, target_genomes)

    def requires(self):
        yield self.clone(PrepareFiles)
        yield self.clone(Chaining)
        yield self.clone(TransMap)
        if self.augustus is True:
            yield self.clone(Augustus)
        #yield EvaluateTranscripts()
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
            yield self.clone(GenomeFasta, **args)
            yield self.clone(GenomeTwoBit, **args)
            yield self.clone(GenomeSizes, **args)
            yield self.clone(GenomeFlatFasta, **args)
            yield self.clone(GenomeFastaIndex, **args)


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
        args = {'annotation_gp': os.path.join(base_dir, annotation + '.gp'),
                 'annotation_attrs': os.path.join(base_dir, annotation + '.attrs.tsv'),
                 'transcript_fasta':  os.path.join(base_dir, annotation + '.fa'),
                 'transcript_flat_fasta': os.path.join(base_dir, annotation + '.fa.flat'),
                 'transcript_bed': os.path.join(base_dir, annotation + '.bed'),
                 'ref_psl': os.path.join(base_dir, annotation + '.psl')}
        return args

    def requires(self):
        args = self.get_args(self.work_dir, self.annotation)
        args.update(GenomeFiles.get_args(self.work_dir, self.ref_genome))
        yield self.clone(Gff3ToGenePred, **args)
        yield self.clone(Gff3ToAttrs, **args)
        yield self.clone(TranscriptBed, **args)
        yield self.clone(TranscriptFasta, **args)
        yield self.clone(FlatTranscriptFasta, **args)
        yield self.clone(FakePsl, **args)


@inherits(ReferenceFiles)
class Gff3ToGenePred(AbstractAtomicFileTask):
    """
    Generates a genePred from a gff3 file.
    """
    annotation_gp = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.annotation_gp)

    def run(self):
        cmd = ['gff3ToGenePred', '-rnaNameAttr=transcript_id', '-geneNameAttr=gene_id', '-honorStartStopCodons',
               self.annotation, '/dev/stdout']
        self.run_cmd(cmd)


@inherits(ReferenceFiles)
class Gff3ToAttrs(luigi.Task):
    """
    Uses the gff3 parser to extract the attributes table.
    """
    annotation_gp = luigi.Parameter()
    annotation_attrs = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.annotation_attrs)

    def run(self):
        results = tools.gff3.extract_attrs(self.annotation)
        with self.output().open('w') as outf:
            results.to_csv(outf, sep='\t', index_label='tx_id')


@requires(Gff3ToGenePred)
class TranscriptBed(AbstractAtomicFileTask):
    """
    Produces a BED record from the input genePred annotation. Makes use of Kent tool genePredToBed
    """
    transcript_bed = luigi.Parameter()
    annotation_gp = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.transcript_bed)

    def run(self):
        cmd = ['genePredToBed', self.annotation_gp, '/dev/stdout']
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
    ref_psl = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.ref_psl)

    def run(self):
        cmd = ['genePredToFakePsl', '-chromSize={}'.format(self.sizes), 'noDB',
               self.annotation_gp, '/dev/stdout', '/dev/null']
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
        job_store = os.path.join(self.work_dir, 'toil', 'augustus', self.genome)
        toil_options = self.prepare_toil_options(job_store)
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
                'transcript_fasta': ref_files['transcript_fasta'], 'ref_psl': ref_files['ref_psl'],
                'annotation_gp': ref_files['annotation_gp'],
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
            yield self.clone(TransMapPsl, tm_args=args, genome=target_genome)
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
        psl_cmd = ['pslMap', '-chainMapFile', self.tm_args['ref_psl'],
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
                outf.write(psl_rec.psl_string() + '\n')


@requires(TransMapPsl)
class TransMapGp(AbstractAtomicFileTask):
    """
    Produces the final transMapped genePred
    """
    def output(self):
        return luigi.LocalTarget(self.tm_args['tm_gp'])

    def run(self):
        cmd = ['transMapPslToGenePred', '-nonCodingGapFillMax=80', '-codingGapFillMax=50',
               self.tm_args['annotation_gp'], self.tm_args['tm_psl'], '/dev/stdout']
        self.run_cmd(cmd)


@inherits(RunCat)
class Augustus(luigi.WrapperTask):
    """
    Runs AugustusTM(R) on the coding output from transMap.
    """
    @staticmethod
    def get_args(work_dir, genome, ref_genome, annotation, augustus_hints_db, tm_cfg, tmr_cfg):
        tm_args = TransMap.get_args(work_dir, genome, ref_genome, annotation)
        tgt_genome_files = GenomeFiles.get_args(work_dir, genome)
        annotation_files = ReferenceFiles.get_args(work_dir, annotation)
        base_dir = os.path.join(work_dir, 'augustus')
        if augustus_hints_db is not None:
            augustus_gp = os.path.join(base_dir, genome + '.TM.gp')
        else:
            augustus_gp = os.path.join(base_dir, genome + '.TMR.gp')
        args = {'ref_genome': ref_genome, 'genome': genome,
                'genome_fasta': os.path.abspath(tgt_genome_files['fasta']),
                'annotation_gp': os.path.abspath(annotation_files['annotation_gp']),
                'ref_psl': os.path.abspath(tm_args['ref_psl']),
                'annotation_attrs': os.path.abspath(annotation_files['annotation_attrs']),
                'tm_gp': os.path.abspath(tm_args['tm_gp']),
                'tm_psl': os.path.abspath(tm_args['tm_psl']),
                'augustus_gp': os.path.abspath(augustus_gp),
                'augustus_hints_db': os.path.abspath(augustus_hints_db) if augustus_hints_db is not None else None,
                'tm_cfg': os.path.abspath(tm_cfg),
                'tmr_cfg': os.path.abspath(tmr_cfg)}
        return args

    def validate(self):
        # TODO: make sure that all input args exist
        pass

    def requires(self):
        self.validate()
        target_genomes = RunCat.resolve_augustus_genomes(self.hal, self.target_genomes, self.augustus_genomes)
        for target_genome in target_genomes:
            args = self.get_args(self.work_dir, target_genome, self.ref_genome, self.annotation, self.augustus_hints_db,
                                 self.tm_cfg, self.tmr_cfg)
            yield self.clone(AugustusDriverTask, augustus_args=args, genome=target_genome)


@inherits(Augustus)
class AugustusDriverTask(tools.toilInterface.ToilTask):
    """
    Task for per-genome launching of a toil pipeline for running Augustus.
    """
    augustus_args = luigi.DictParameter()  # dict to pass directory to TMR toil module
    genome = luigi.Parameter()
    
    def output(self):
        return luigi.LocalTarget(self.augustus_args['augustus_gp'])
    
    def requires(self):
        return self.clone(TransMap)

    def extract_coding_genes(self):
        """extracts only coding genes from the input genePred, returning a path to a tmp file"""
        coding_gp = tools.fileOps.get_tmp_file()
        attrs = tools.misc.read_attributes_tsv(self.augustus_args['annotation_attrs'])
        names = set(attrs[attrs.tx_biotype == 'protein_coding'].index)
        with open(coding_gp, 'w') as outf:
            for name, tx in tools.transcripts.gene_pred_iterator(self.augustus_args['tm_gp']):
                if tools.psl.strip_alignment_numbers(name) in names:
                    outf.write(tx.get_gene_pred() + '\n')
        if os.path.getsize(coding_gp) == 0:
            raise RuntimeError('Unable to extract coding transcripts from the transMap genePred.')
        return coding_gp

    def convert_gtf_gp(self, gtf_path):
        """converts the Augustus output GTF to genePred"""
        with self.output().open('w') as outf:
            cmd = ['gtfToGenePred', gtf_path, '/dev/stdout']
            tools.procOps.run_proc(cmd, stdout=outf)

    def run(self):
        job_store = os.path.join(self.work_dir, 'toil', 'augustus', self.genome)
        toil_options = self.prepare_toil_options(job_store)
        coding_gp = self.extract_coding_genes()
        augustus_results = augustus(self.augustus_args, coding_gp, toil_options)
        os.remove(coding_gp)
        with tools.fileOps.TemporaryFilePath() as tmp_gtf:
            with open(tmp_gtf, 'w') as outf:
                tools.fileOps.print_rows(outf, augustus_results)
            self.convert_gtf_gp(tmp_gtf)


class EvaluateTranscripts(luigi.WrapperTask):
    """
    Evaluates all transcripts for important features.
    """
    def validate(self):
        # TODO: make sure that all input args exist
        pass

    @staticmethod
    def get_args(work_dir, genome, ref_genome, annotation, augustus_hints_db, tm_cfg, tmr_cfg, augustus):
        tm_args = TransMap.get_args(work_dir, genome, ref_genome, annotation)
        tgt_genome_files = GenomeFiles.get_args(work_dir, genome)
        annotation_files = ReferenceFiles.get_args(work_dir, annotation)
        base_dir = os.path.join(work_dir, 'tm_eval')
        args = {'db': os.path.join(base_dir, 'evaluations.db'),
                'tm_psl': tm_args['tm_psl'],
                'tm_gp': tm_args['tm_gp'],
                'annotation_gp': annotation_files['annotation_gp'],
                'annotation_attrs': annotation_files['annotation_attrs'],
                'genome_fasta': tgt_genome_files['fasta'],
                'genome': genome,
                'classify_table': genome + '_Classify',
                'details_table': genome + '_Details'}
        if augustus is True:
            aug_args = Augustus.get_args(work_dir, genome, ref_genome, annotation, augustus_hints_db, tm_cfg, tmr_cfg)
            args['augustus_gp'] = aug_args['augustus_gp']
        else:
            args['augustus_gp'] = None
        return args

    def requires(self):
        self.validate()
        target_genomes = RunCat.resolve_augustus_genomes(self.hal, self.target_genomes, self.augustus_genomes)
        for target_genome in target_genomes:
            args = self.get_args(self.work_dir, target_genome, self.ref_genome, self.annotation, self.augustus_hints_db,
                                 self.tm_cfg, self.tmr_cfg, self.augustus)
            yield self.clone(EvaluateDriverTask, eval_args=args, genome=target_genome)


@inherits(EvaluateTranscripts)
class EvaluateDriverTask(tools.toilInterface.ToilTask):
    """
    Task for per-genome launching of a toil pipeline for running Augustus.
    """
    eval_args = luigi.DictParameter()  # dict to pass directory to evaluate toil module
    genome = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.augustus_args['augustus_gp'])

    def requires(self):
        return self.clone(TransMap)

    def convert_gtf_gp(self, gtf_path):
        """converts the Augustus output GTF to genePred"""
        with self.output().open('w') as outf:
            cmd = ['gtfToGenePred', gtf_path, '/dev/stdout']
            tools.procOps.run_proc(cmd, stdout=outf)

    def run(self):
        job_store = os.path.join(self.work_dir, 'toil', 'augustus', self.genome)
        toil_options = self.prepare_toil_options(job_store)
        coding_gp = self.extract_coding_genes()
        augustus_results = augustus(self.augustus_args, coding_gp, toil_options)
        os.remove(coding_gp)
        with tools.fileOps.TemporaryFilePath() as tmp_gtf:
            with open(tmp_gtf, 'w') as outf:
                tools.fileOps.print_rows(outf, augustus_results)
            self.convert_gtf_gp(tmp_gtf)


if __name__ == '__main__':
    luigi.build([Augustus(hal='1509.hal',
                target_genomes=('C57BL_6NJ',), work_dir='cat_work', ref_genome='C57B6J', augustus=True,
                annotation='Mus_musculus.GRCm38.83.gff3')], logging_conf_file='logging.cfg')
