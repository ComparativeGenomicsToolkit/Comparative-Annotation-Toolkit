#!/usr/bin/env python2.7
"""
Comparative Annotation Toolkit.
"""
import argparse
import itertools
import logging
import multiprocessing
import os
import tempfile
from collections import OrderedDict

import luigi
import luigi.contrib.sqla
from luigi.util import requires

import tools.bio
import tools.fileOps
import tools.gff3
import tools.hal
import tools.misc
import tools.nameConversions
import tools.procOps
import tools.mathOps
import tools.psl
import tools.sqlInterface
import tools.sqlite
import tools.hintsDatabaseInterface
import tools.transcripts
from tools.luigiAddons import multiple_requires
from align_transcripts import align_transcripts
from augustus import augustus
from augustus_cgp import augustus_cgp
from base_tasks import PipelineTask, PipelineWrapperTask, ToilTask, AbstractAtomicFileTask
from chaining import chaining
from classify import classify
from consensus import generate_consensus
from filter_transmap import filter_transmap
from hgm import hgm, parse_hgm_gtf
from transmap_classify import transmap_classify
from plots import generate_plots

logger = logging.getLogger(__name__)


###
# Pipeline exceptions
###


class UserException(Exception):
    """generic exception to use when a user makes a mistake"""
    pass


class ToolMissingException(UserException):
    """exception to use when a tool is missing, usually checked in a task validate() method"""
    pass


class InputMissingException(UserException):
    """exception to use when input data are missing"""
    pass


class InvalidInputException(UserException):
    """exception to use when something about the input is invalid"""
    pass


###
# pipeline tasks
###


class RunCat(PipelineWrapperTask):
    """
    Task that executes the entire pipeline.
    """
    def validate(self):
        """General input validation"""
        pipeline_args = self.get_pipeline_args()
        if not os.path.exists(pipeline_args.hal):
            raise InputMissingException('HAL file not found at {}.'.format(pipeline_args.hal))
        for d in [pipeline_args.out_dir, pipeline_args.work_dir]:
            if not os.path.exists(d):
                if not tools.fileOps.dir_is_writeable(os.path.dirname(d)):
                    raise UserException('Cannot create directory {}.'.format(d))
            else:
                if not tools.fileOps.dir_is_writeable(d):
                    raise UserException('Directory {} is not writeable.'.format(d))
        if not os.path.exists(pipeline_args.annotation):
            raise InputMissingException('Annotation file {} not found.'.format(pipeline_args.annotation))
        if pipeline_args.augustus_hints_db is not None and not os.path.exists(pipeline_args.augustus_hints_db):
            raise InputMissingException('Augustus hints DB {} not found.'.format(pipeline_args.augustus_hints_db))
        # TODO: validate augustus species, tm/tmr/cgp/param files.
        if pipeline_args.ref_genome not in pipeline_args.hal_genomes:
            raise InvalidInputException('Reference genome {} not present in HAL.'.format(pipeline_args.ref_genome))
        missing_genomes = {g for g in pipeline_args.target_genomes if g not in pipeline_args.hal_genomes}
        if len(missing_genomes) > 0:
            missing_genomes = ','.join(missing_genomes)
            raise InvalidInputException('Target genomes {} not present in HAL.'.format(missing_genomes))
        if pipeline_args.ref_genome in pipeline_args.target_genomes:
            raise InvalidInputException('A target genome cannot be the reference genome.')
        return pipeline_args

    def requires(self):
        pipeline_args = self.validate()
        yield self.clone(PrepareFiles)
        yield self.clone(Chaining)
        yield self.clone(TransMap)
        yield self.clone(EvaluateTransMap)
        yield self.clone(FilterTransMap)
        if self.augustus is True:
            yield self.clone(Augustus)
        if self.augustus_cgp is True:
            yield self.clone(AugustusCgp)
        if pipeline_args.augustus_tmr is True:
            yield self.clone(Hgm)
        yield self.clone(AlignTranscripts)
        yield self.clone(EvaluateTranscripts)
        yield self.clone(Consensus)
        yield self.clone(Plots)


class PrepareFiles(PipelineWrapperTask):
    """
    Wrapper for file preparation tasks GenomeFiles and ReferenceFiles
    """
    def requires(self):
        yield self.clone(GenomeFiles)
        yield self.clone(ReferenceFiles)


class GenomeFiles(PipelineWrapperTask):
    """
    WrapperTask for producing all genome files.

    GenomeFiles -> GenomeFasta -> GenomeTwoBit -> GenomeFlatFasta -> GenomeFastaIndex
                -> GenomeSizes

    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'genome_files')
        args = argparse.Namespace()
        args.genome = genome
        args.fasta = os.path.join(base_dir, genome + '.fa')
        args.two_bit = os.path.join(base_dir, genome + '.2bit')
        args.sizes = os.path.join(base_dir, genome + '.chrom.sizes')
        args.flat_fasta = os.path.join(base_dir, genome + '.fa.flat')
        return args

    def validate(self):
        for haltool in ['hal2fasta', 'halStats']:
            if not tools.misc.is_exec(haltool):
                    raise ToolMissingException('{} from the HAL tools package not in global path'.format(haltool))
        if not tools.misc.is_exec('faToTwoBit'):
            raise ToolMissingException('faToTwoBit tool from the Kent tools package not in global path.')
        if not tools.misc.is_exec('pyfasta'):
            raise ToolMissingException('pyfasta wrapper not found in global path.')

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for genome in itertools.chain(pipeline_args.target_genomes, [self.ref_genome]):
            args = self.get_args(pipeline_args, genome)
            yield self.clone(GenomeFasta, **vars(args))
            yield self.clone(GenomeTwoBit, **vars(args))
            yield self.clone(GenomeSizes, **vars(args))
            yield self.clone(GenomeFlatFasta, **vars(args))


class GenomeFasta(AbstractAtomicFileTask):
    """
    Produce a fasta file from a hal file. Requires hal2fasta.
    """
    genome = luigi.Parameter()
    fasta = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.fasta)

    def run(self):
        logger.info('Extracting fasta for {}.'.format(self.genome))
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
        logger.info('Converting fasta for {} to 2bit.'.format(self.genome))
        cmd = ['faToTwoBit', self.fasta, '/dev/stdout']
        self.run_cmd(cmd)


class GenomeSizes(AbstractAtomicFileTask):
    """
    Produces a genome chromosome sizes file. Requires halStats.
    """
    genome = luigi.Parameter()
    sizes = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.sizes)

    def run(self):
        logger.info('Extracting chromosome sizes for {}.'.format(self.genome))
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
        logger.info('Flattening fasta for {}.'.format(self.genome))
        cmd = ['pyfasta', 'flatten', self.fasta]
        tools.procOps.run_proc(cmd)


class ReferenceFiles(PipelineWrapperTask):
    """
    WrapperTask for producing annotation files.

    ReferenceFiles -> Gff3ToGenePred -> TranscriptBed -> TranscriptFasta -> FlatTranscriptFasta
                            V
                         FakePsl
    """
    @staticmethod
    def get_args(pipeline_args):
        base_dir = os.path.join(pipeline_args.work_dir, 'reference')
        annotation = os.path.splitext(os.path.basename(pipeline_args.annotation))[0]
        args = argparse.Namespace()
        args.annotation_gp = os.path.join(base_dir, annotation + '.gp')
        args.transcript_fasta = os.path.join(base_dir, annotation + '.fa')
        args.transcript_flat_fasta = os.path.join(base_dir, annotation + '.fa.flat')
        args.transcript_bed = os.path.join(base_dir, annotation + '.bed')
        args.ref_psl = os.path.join(base_dir, annotation + '.psl')
        args.__dict__.update(**vars(GenomeFiles.get_args(pipeline_args, pipeline_args.ref_genome)))
        return args

    def validate(self):
        for tool in ['gff3ToGenePred', 'genePredToBed', 'genePredToFakePsl']:
            if not tools.misc.is_exec(tool):
                    raise ToolMissingException('{} from the Kent tools package not in global path'.format(tool))

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        args = self.get_args(pipeline_args)
        yield self.clone(Gff3ToGenePred, **vars(args))
        yield self.clone(Gff3ToAttrs, **vars(args))
        yield self.clone(TranscriptBed, **vars(args))
        yield self.clone(TranscriptFasta, **vars(args))
        yield self.clone(FlatTranscriptFasta, **vars(args))
        yield self.clone(FakePsl, **vars(args))


class Gff3ToGenePred(AbstractAtomicFileTask):
    """
    Generates a genePred from a gff3 file.
    """
    annotation_gp = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.annotation_gp)

    def run(self):
        logger.info('Converting annotation gff3 to genePred.')
        cmd = ['gff3ToGenePred', '-rnaNameAttr=transcript_id', '-geneNameAttr=gene_id', '-honorStartStopCodons',
               self.annotation, '/dev/stdout']
        self.run_cmd(cmd)


class Gff3ToAttrs(PipelineTask):
    """
    Uses the gff3 parser to extract the attributes table, converting the table into a sqlite database.
    """
    table = tools.sqlInterface.Annotation.__tablename__

    def output(self):
        pipeline_args = self.get_pipeline_args()
        database = self.__class__.get_database(pipeline_args, pipeline_args.ref_genome)
        tools.fileOps.ensure_file_dir(database)
        conn_str = 'sqlite:///{}'.format(database)
        digest = tools.fileOps.hashfile(self.annotation)
        attrs_table = luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                          target_table=self.table,
                                                          update_id='_'.join([self.table, digest]))
        return attrs_table

    def run(self):
        logger.info('Extracting gff3 attributes to sqlite database.')
        results = tools.gff3.extract_attrs(self.annotation)
        if 'protein_coding' not in results.TranscriptBiotype[1] or 'protein_coding' not in results.GeneBiotype[1]:
            logger.warning('No protein_coding annotations found!')
        pipeline_args = self.get_pipeline_args()
        database = self.__class__.get_database(pipeline_args, pipeline_args.ref_genome)
        with tools.sqlite.ExclusiveSqlConnection(database) as engine:
            results.to_sql(self.table, engine, if_exists='replace')
        self.output().touch()


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
        logger.info('Converting annotation genePred to BED.')
        cmd = ['genePredToBed', self.annotation_gp, '/dev/stdout']
        self.run_cmd(cmd)


@multiple_requires(GenomeFlatFasta, TranscriptBed)
class TranscriptFasta(AbstractAtomicFileTask):
    """
    Produces a fasta for each transcript.
    """
    transcript_fasta = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.transcript_fasta)

    def run(self):
        logger.info('Extracting reference annotation fasta.')
        seq_dict = tools.bio.get_sequence_dict(self.fasta, upper=False)
        seqs = {tx.name: tx.get_mrna(seq_dict) for tx in tools.transcripts.transcript_iterator(self.transcript_bed)}
        with self.output().open('w') as outf:
            for name, seq in seqs.iteritems():
                tools.bio.write_fasta(outf, name, seq)


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
        logger.info('Flattening reference annotation fasta.')
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
        logger.info('Generating annotation fake PSL.')
        cmd = ['genePredToFakePsl', '-chromSize={}'.format(self.sizes), 'noDB',
               self.annotation_gp, '/dev/stdout', '/dev/null']
        self.run_cmd(cmd)


class Chaining(ToilTask):
    """
    Task that launches the Chaining toil pipeline. This pipeline operates on all genomes at once to reduce the
    repeated downloading of the HAL file.
    """

    @staticmethod
    def get_args(pipeline_args):
        base_dir = os.path.join(pipeline_args.work_dir, 'chaining')
        ref_files = GenomeFiles.get_args(pipeline_args, pipeline_args.ref_genome)
        tgt_files = {genome: GenomeFiles.get_args(pipeline_args, genome) for genome in pipeline_args.target_genomes}
        tgt_two_bits = {genome: tgt_files[genome].two_bit for genome in pipeline_args.target_genomes}
        chain_files = {genome: os.path.join(base_dir, '{}-{}.chain'.format(pipeline_args.ref_genome, genome))
                       for genome in pipeline_args.target_genomes}
        args = argparse.Namespace()
        args.hal = pipeline_args.hal
        args.ref_genome = pipeline_args.ref_genome
        args.query_two_bit = ref_files.two_bit
        args.query_sizes = ref_files.sizes
        args.target_two_bits = tgt_two_bits
        args.chain_files = chain_files
        return args

    def output(self):
        pipeline_args = self.get_pipeline_args()
        chain_args = self.get_args(pipeline_args)
        for path in chain_args.chain_files.itervalues():
            yield luigi.LocalTarget(path)

    def validate(self):
        if not tools.misc.is_exec('halLiftover'):
            raise ToolMissingException('halLiftover from the halTools package not in global path.')
        for tool in ['pslPosTarget', 'axtChain', 'chainMergeSort']:
            if not tools.misc.is_exec(tool):
                    raise ToolMissingException('{} from the Kent tools package not in global path.'.format(tool))

    def requires(self):
        yield self.clone(PrepareFiles)

    def run(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        logger.info('Launching Pairwise Chaining toil pipeline.')
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'chaining')
        toil_options = self.prepare_toil_options(toil_work_dir)
        chain_args = self.get_args(pipeline_args)
        chaining(chain_args, toil_options)
        logger.info('Pairwise Chaining toil pipeline is complete.')


class TransMap(PipelineWrapperTask):
    """
    Runs transMap.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'transMap')
        ref_files = ReferenceFiles.get_args(pipeline_args)
        tgt_genome_files = GenomeFiles.get_args(pipeline_args, genome)
        chain_args = Chaining.get_args(pipeline_args)
        args = argparse.Namespace()
        args.two_bit = tgt_genome_files.two_bit
        args.chain_file = chain_args.chain_files[genome]
        args.transcript_fasta = ref_files.transcript_fasta
        args.ref_psl = ref_files.ref_psl
        args.annotation_gp = ref_files.annotation_gp
        args.tm_psl = os.path.join(base_dir, genome + '.psl')
        args.tm_gp = os.path.join(base_dir, genome + '.gp')
        args.tm_gtf = os.path.join(base_dir, genome + '.gtf')
        return args

    def validate(self):
        for tool in ['pslMap', 'pslRecalcMatch', 'transMapPslToGenePred']:
            if not tools.misc.is_exec(tool):
                    raise ToolMissingException('{} from the Kent tools package not in global path.'.format(tool))

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            yield self.clone(TransMapPsl, genome=target_genome)
            yield self.clone(TransMapGp, genome=target_genome)
            yield self.clone(TransMapGtf, genome=target_genome)


class TransMapPsl(PipelineTask):
    """
    Runs transMap. Requires Kent tools pslMap, postTransMapChain, pslRecalcMatch
    """
    genome = luigi.Parameter()

    def output(self):
        tm_args = self.get_module_args(TransMap, genome=self.genome)
        return luigi.LocalTarget(tm_args.tm_psl)

    def requires(self):
        return self.clone(PrepareFiles), self.clone(Chaining)

    def run(self):
        logger.info('Running transMap for {}.'.format(self.genome))
        tm_args = self.get_module_args(TransMap, genome=self.genome)
        cmd = [['pslMap', '-chainMapFile', tm_args.ref_psl, tm_args.chain_file, '/dev/stdout'],
               ['postTransMapChain', '/dev/stdin', '/dev/stdout'],
               ['sort', '-k', '14,14', '-k', '16,16n'],
               ['pslRecalcMatch', '/dev/stdin', tm_args.two_bit, tm_args.transcript_fasta, 'stdout'],
               ['pslCDnaFilter', '-localNearBest=0.0001', '-minCover=0.1', '/dev/stdin', '/dev/stdout'],
               ['awk', '$17 - $16 < 3000000 {print $0}']]  # hard coded filter for 3mb transcripts
        # hacky way to make unique - capture output to a file, then process
        tmp_file = luigi.LocalTarget(is_tmp=True)
        with tmp_file.open('w') as tmp_fh:
            tools.procOps.run_proc(cmd, stdout=tmp_fh, stderr='/dev/null')
        tools.fileOps.ensure_file_dir(self.output().path)
        with self.output().open('w') as outf:
            for psl_rec in tools.psl.psl_iterator(tmp_file.path, make_unique=True):
                outf.write('\t'.join(psl_rec.psl_string()) + '\n')


@requires(TransMapPsl)
class TransMapGp(AbstractAtomicFileTask):
    """
    Produces the final transMapped genePred
    """
    def output(self):
        tm_args = self.get_module_args(TransMap, genome=self.genome)
        return luigi.LocalTarget(tm_args.tm_gp)

    def run(self):
        tm_args = self.get_module_args(TransMap, genome=self.genome)
        logger.info('Converting transMap PSL to genePred for {}.'.format(self.genome))
        cmd = ['transMapPslToGenePred', '-nonCodingGapFillMax=80', '-codingGapFillMax=50',
               tm_args.annotation_gp, tm_args.tm_psl, '/dev/stdout']
        self.run_cmd(cmd)


@requires(TransMapGp)
class TransMapGtf(PipelineTask):
    """
    Converts the transMap genePred to GTF
    """
    def output(self):
        tm_args = self.get_module_args(TransMap, genome=self.genome)
        return luigi.LocalTarget(tm_args.tm_gtf)

    def run(self):
        tm_args = self.get_module_args(TransMap, genome=self.genome)
        logger.info('Converting transMap genePred to GTF for {}.'.format(self.genome))
        tools.misc.convert_gp_gtf(self.output(), luigi.LocalTarget(tm_args.tm_gp))


class EvaluateTransMap(PipelineWrapperTask):
    """
    Evaluates transMap alignments.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        tm_args = TransMap.get_args(pipeline_args, genome)
        args = argparse.Namespace()
        args.db_path = pipeline_args.dbs[genome]
        args.tm_psl = tm_args.tm_psl
        args.ref_psl = ReferenceFiles.get_args(pipeline_args).ref_psl
        args.tm_gp = tm_args.tm_gp
        args.annotation_gp = tm_args.annotation_gp
        args.annotation_gp = ReferenceFiles.get_args(pipeline_args).annotation_gp
        args.genome = genome
        args.fasta = GenomeFiles.get_args(pipeline_args, genome).fasta
        args.ref_genome = pipeline_args.ref_genome
        return args

    def validate(self):
        pass

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            tm_eval_args = EvaluateTransMap.get_args(pipeline_args, target_genome)
            yield self.clone(EvaluateTransMapDriverTask, tm_eval_args=tm_eval_args, genome=target_genome)


class EvaluateTransMapDriverTask(PipelineTask):
    """
    Task for per-genome launching of a toil pipeline for aligning transcripts to their parent.
    """
    genome = luigi.Parameter()
    tm_eval_args = luigi.Parameter()
    table = tools.sqlInterface.TmEval.__tablename__

    def write_to_sql(self, df):
        """Load the results into the SQLite database"""
        with tools.sqlite.ExclusiveSqlConnection(self.tm_eval_args.db_path) as engine:
            df.to_sql(self.table, engine, if_exists='replace')
            self.output().touch()
            logger.info('Loaded table: {}.{}'.format(self.genome, self.table))

    def output(self):
        pipeline_args = self.get_pipeline_args()
        tools.fileOps.ensure_file_dir(self.tm_eval_args.db_path)
        conn_str = 'sqlite:///{}'.format(self.tm_eval_args.db_path)
        return luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                   target_table=self.table,
                                                   update_id='_'.join([self.table, str(hash(pipeline_args))]))

    def requires(self):
        if self.no_evaluate_dependency is True:
            return
        return self.clone(TransMap), self.clone(ReferenceFiles)

    def run(self):
        logger.info('Evaluating transMap results for {}.'.format(self.genome))
        results = transmap_classify(self.tm_eval_args)
        self.write_to_sql(results)


class FilterTransMap(PipelineWrapperTask):
    """
    Filters transMap alignments for paralogs, as well as multiple chromosomes if the --resolve-split-genes flag is set.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'filtered_transMap')
        args = argparse.Namespace()
        args.tm_gp = TransMap.get_args(pipeline_args, genome).tm_gp
        args.filtered_tm_gp = os.path.join(base_dir, genome + '.filtered.gp')
        args.filtered_tm_gtf = os.path.join(base_dir, genome + '.filtered.gtf')
        args.db_path = pipeline_args.dbs[genome]
        args.ref_db_path = pipeline_args.dbs[pipeline_args.ref_genome]
        args.resolve_split_genes = pipeline_args.resolve_split_genes
        args.metrics_json = os.path.join(PipelineTask.get_metrics_dir(pipeline_args, genome), 'filter_tm_metrics.json')
        return args

    def validate(self):
        pass

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            filter_tm_args = FilterTransMap.get_args(pipeline_args, target_genome)
            yield self.clone(FilterTransMapDriverTask, filter_tm_args=filter_tm_args, genome=target_genome)


class FilterTransMapDriverTask(PipelineTask):
    """
    Driver task for per-genome transMap filtering.
    """
    genome = luigi.Parameter()
    filter_tm_args = luigi.Parameter()
    eval_table = tools.sqlInterface.TmFilterEval.__tablename__
    cutoff_table = tools.sqlInterface.TmFit.__tablename__

    def write_to_sql(self, updated_df, fit_df, filter_table_target, fit_table_target):
        """Load the results into the SQLite database"""
        with tools.sqlite.ExclusiveSqlConnection(self.filter_tm_args.db_path) as engine:
            updated_df.to_sql(self.eval_table, engine, if_exists='replace')
            filter_table_target.touch()
            logger.info('Loaded table: {}.{}'.format(self.genome, self.eval_table))
            fit_df.to_sql(self.cutoff_table, engine, if_exists='replace')
            fit_table_target.touch()
            logger.info('Loaded table: {}.{}'.format(self.genome, self.cutoff_table))

    def output(self):
        pipeline_args = self.get_pipeline_args()
        tools.fileOps.ensure_file_dir(self.filter_tm_args.db_path)
        conn_str = 'sqlite:///{}'.format(self.filter_tm_args.db_path)
        return (luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                    target_table=self.eval_table,
                                                    update_id='_'.join([self.eval_table, str(hash(pipeline_args))])),
                luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                    target_table=self.cutoff_table,
                                                    update_id='_'.join([self.cutoff_table, str(hash(pipeline_args))])),
                luigi.LocalTarget(self.filter_tm_args.filtered_tm_gp),
                luigi.LocalTarget(self.filter_tm_args.filtered_tm_gtf),
                luigi.LocalTarget(self.filter_tm_args.metrics_json))

    def requires(self):
        if self.no_evaluate_dependency is True:
            return
        return self.clone(EvaluateTransMap)

    def run(self):
        logger.info('Filtering transMap results for {}.'.format(self.genome))
        filter_table_target, fit_table_target, filtered_tm_gp, filtered_tm_gtf, metrics_json = self.output()
        metrics_dict, updated_df, fit_df = filter_transmap(self.filter_tm_args, filtered_tm_gp)
        PipelineTask.write_metrics(metrics_dict, metrics_json)
        tools.misc.convert_gp_gtf(filtered_tm_gtf, filtered_tm_gp)
        self.write_to_sql(updated_df, fit_df, filter_table_target, fit_table_target)


class Augustus(PipelineWrapperTask):
    """
    Runs AugustusTM(R) on the coding output from transMap.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'augustus')
        args = argparse.Namespace()
        args.ref_genome = pipeline_args.ref_genome
        args.genome = genome
        args.genome_fasta = GenomeFiles.get_args(pipeline_args, genome).fasta
        args.annotation_gp = ReferenceFiles.get_args(pipeline_args).annotation_gp
        args.ref_db_path = PipelineTask.get_database(pipeline_args, pipeline_args.ref_genome)
        args.filtered_tm_gp = FilterTransMap.get_args(pipeline_args, genome).filtered_tm_gp
        tm_args = TransMap.get_args(pipeline_args, genome)
        args.ref_psl = tm_args.ref_psl
        args.tm_psl = tm_args.tm_psl
        args.augustus_tm_gp = os.path.join(base_dir, genome + '.augTM.gp')
        args.augustus_tm_gtf = os.path.join(base_dir, genome + '.augTM.gtf')
        args.augustus_hints_db = pipeline_args.augustus_hints_db if pipeline_args.augustus_tmr else None
        args.tm_cfg = pipeline_args.tm_cfg
        args.tmr_cfg = pipeline_args.tmr_cfg
        args.augustus_species = pipeline_args.augustus_species
        if pipeline_args.augustus_tmr:
            args.augustus_tmr = True
            args.augustus_tmr_gp = os.path.join(base_dir, genome + '.augTMR.gp')
            args.augustus_tmr_gtf = os.path.join(base_dir, genome + '.augTMR.gtf')
        return args

    def validate(self):
        if not tools.misc.is_exec('augustus'):
            raise ToolMissingException('auxiliary program augustus not in global path.')
        if not tools.misc.is_exec('transMap2hints.pl'):
            raise ToolMissingException('auxiliary program transMap2hints.pl from the Augustus package '
                                       'not in global path.')

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            yield self.clone(AugustusDriverTask, genome=target_genome)


class AugustusDriverTask(ToilTask):
    """
    Task for per-genome launching of a toil pipeline for running Augustus.
    """
    genome = luigi.Parameter()

    def output(self):
        pipeline_args = self.get_pipeline_args()
        augustus_args = Augustus.get_args(pipeline_args, self.genome)
        yield luigi.LocalTarget(augustus_args.augustus_tm_gp)
        yield luigi.LocalTarget(augustus_args.augustus_tm_gtf)
        if augustus_args.augustus_tmr:
            yield luigi.LocalTarget(augustus_args.augustus_tmr_gp)
            yield luigi.LocalTarget(augustus_args.augustus_tmr_gtf)

    def requires(self):
        if self.no_evaluate_dependency is True:
            return
        return self.clone(FilterTransMap)

    def extract_coding_genes(self, augustus_args):
        """extracts only coding genes from the input genePred, returning a path to a tmp file"""
        coding_gp = tools.fileOps.get_tmp_file()
        attrs = tools.sqlInterface.read_attrs(augustus_args.ref_db_path)
        names = set(attrs[attrs.TranscriptBiotype == 'protein_coding'].index)
        with open(coding_gp, 'w') as outf:
            for tx in tools.transcripts.gene_pred_iterator(augustus_args.filtered_tm_gp):
                if tools.nameConversions.strip_alignment_numbers(tx.name) in names:
                    tools.fileOps.print_row(outf, tx.get_gene_pred())
        if os.path.getsize(coding_gp) == 0:
            raise InvalidInputException('Unable to extract coding transcripts from the filtered transMap genePred.')
        return coding_gp

    def run(self):
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'augustus', self.genome)
        logger.info('Launching Augustus toil pipeline on {}.'.format(self.genome))
        toil_options = self.prepare_toil_options(toil_work_dir)
        augustus_args = self.get_module_args(Augustus, genome=self.genome)
        coding_gp = self.extract_coding_genes(augustus_args)
        augustus(augustus_args, coding_gp, toil_options)
        logger.info('Augustus toil pipeline for {} completed.'.format(self.genome))
        os.remove(coding_gp)
        for out_gp, out_gtf in tools.misc.pairwise(self.output()):
            tools.misc.convert_gtf_gp(out_gp, out_gtf)


class AugustusCgp(ToilTask):
    """
    Task for launching the AugustusCGP toil pipeline
    """
    @staticmethod
    def get_args(pipeline_args):
        # add reference to the target genomes
        tgt_genomes = list(pipeline_args.target_genomes) + [pipeline_args.ref_genome]
        fasta_files = {genome: GenomeFiles.get_args(pipeline_args, genome).fasta for genome in tgt_genomes}
        base_dir = os.path.join(pipeline_args.work_dir, 'augustus_cgp')
        # output
        output_gp_files = {genome: os.path.join(base_dir, genome + '.augCGP.gp') for genome in tgt_genomes}
        output_gtf_files = {genome: os.path.join(base_dir, genome + '.augCGP.gtf') for genome in tgt_genomes}
        # transMap files used for assigning parental gene
        tm_gp_files = {genome: FilterTransMap.get_args(pipeline_args, genome).filtered_tm_gp
                       for genome in pipeline_args.target_genomes}
        # add the reference annotation as a pseudo-transMap to assign parents in reference
        tm_gp_files[pipeline_args.ref_genome] = ReferenceFiles.get_args(pipeline_args).annotation_gp
        hints_db = pipeline_args.augustus_hints_db if pipeline_args.augustus_hints_db else None
        args = argparse.Namespace()
        args.genomes = tgt_genomes
        args.fasta_files = fasta_files
        args.tm_gps = tm_gp_files
        args.hal = pipeline_args.hal
        args.ref_genome = pipeline_args.ref_genome
        args.augustus_cgp_gp = output_gp_files
        args.augustus_cgp_gtf = output_gtf_files
        args.species = pipeline_args.augustus_species
        args.chunksize = pipeline_args.maf_chunksize
        args.overlap = pipeline_args.maf_overlap
        args.cgp_param = pipeline_args.cgp_param
        args.hints_db = hints_db
        args.ref_db_path = PipelineTask.get_database(pipeline_args, pipeline_args.ref_genome)
        args.query_sizes = GenomeFiles.get_args(pipeline_args, pipeline_args.ref_genome).sizes
        return args

    def output(self):
        pipeline_args = self.get_pipeline_args()
        cgp_args = self.get_args(pipeline_args)
        for path_dict in [cgp_args.augustus_cgp_gp, cgp_args.augustus_cgp_gtf]:
            for path in path_dict.itervalues():
                yield luigi.LocalTarget(path)

    def validate(self):
        if not tools.misc.is_exec('joingenes'):
            raise ToolMissingException('auxiliary program joingenes from the Augustus package not in global path.')
        if not tools.misc.is_exec('augustus'):
            raise ToolMissingException('augustus not in global path.')
        if not tools.misc.is_exec('hal2maf'):
            raise ToolMissingException('hal2maf from the halTools package not in global path.')
        if not tools.misc.is_exec('gtfToGenePred'):
            raise ToolMissingException('gtfToGenePred from the Kent package not in global path.')
        if not tools.misc.is_exec('genePredToGtf'):
            raise ToolMissingException('genePredToGtf from the Kent package not in global path.')
        if not tools.misc.is_exec('bedtools'):
            raise ToolMissingException('bedtools not in global path.')

    def requires(self):
        self.validate()
        if self.no_evaluate_dependency is True:
            return
        yield self.clone(FilterTransMap), self.clone(Gff3ToAttrs)

    def prepare_cfg(self, pipeline_args):
        """use the config template to create a config file"""
        template = open(pipeline_args.augustus_cgp_cfg_template).read()
        cfg = template.format(ref_genome=pipeline_args.ref_genome,
                              target_genomes=' '.join(pipeline_args.target_genomes))
        out_path = tools.fileOps.get_tmp_file()
        with open(out_path, 'w') as outf:
            outf.write(cfg)
        return out_path

    def evaluate_database(self, pipeline_args):
        """warn the user about the database not containing things"""
        if tools.hintsDatabaseInterface.hints_db_has_annotation(pipeline_args.augustus_hints_db,
                                                                pipeline_args.ref_genome) is False:
            logger.warning('AugustusCGP is being ran without any annotation hints!')
        else:
            logger.info('Hints database has annotation hints for {}.'.format(pipeline_args.ref_genome))
        genomes_with_hints = []
        for genome in pipeline_args.target_genomes:
            if tools.hintsDatabaseInterface.hints_db_has_rnaseq(pipeline_args.augustus_hints_db, genome) is True:
                genomes_with_hints.append(genome)
        if len(genomes_with_hints) == 0:
            logger.warning('AugustusCGP is being ran without any RNA-seq hints!')
        else:
            logger.info('RNA-seq hints found for the following genomes: {}.'.format(','.join(genomes_with_hints)))

    def run(self):
        pipeline_args = self.get_pipeline_args()
        if pipeline_args.augustus_hints_db is None:
            # TODO: if no database is given (de novo setting),
            # create a new DB with the flattened genomes from the HAL alignment
            raise UserException('Cannot run AugustusCGP without a hints database.')
        self.evaluate_database(pipeline_args)
        logger.info('Launching AugustusCGP toil pipeline.')
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'augustus_cgp')
        toil_options = self.prepare_toil_options(toil_work_dir)
        cgp_args = self.get_args(pipeline_args)
        cgp_args.cgp_cfg = self.prepare_cfg(pipeline_args)
        augustus_cgp(cgp_args, toil_options)
        logger.info('AugustusCGP toil pipeline completed.')
        # convert each to genePred as well
        for genome in itertools.chain(pipeline_args.target_genomes, [pipeline_args.ref_genome]):
            gp_target = luigi.LocalTarget(cgp_args.augustus_cgp_gp[genome])
            gtf_target = luigi.LocalTarget(cgp_args.augustus_cgp_gtf[genome])
            tools.misc.convert_gtf_gp(gp_target, gtf_target)
        logger.info('Finished converting AugustusCGP output.')


class Hgm(PipelineWrapperTask):
    """
    Task for launching the HomGeneMapping toil pipeline. This pipeline finds the intron support vector across all
    species in the alignment with RNAseq in the database. It will be launched once for each of transMap, AugustusTM,
    AugustusTMR, AugustusCGP
    """
    @staticmethod
    def get_args(pipeline_args, mode):
        base_dir = os.path.join(pipeline_args.work_dir, 'hgm', mode)
        if mode == 'augCGP':
            # add reference to the target genomes
            tgt_genomes = list(pipeline_args.target_genomes) + [pipeline_args.ref_genome]
            gtf_in_files = {genome: AugustusCgp.get_args(pipeline_args).augustus_cgp_gtf[genome]
                            for genome in tgt_genomes}
        elif mode == 'augTM':
            tgt_genomes = pipeline_args.target_genomes
            gtf_in_files = {genome: Augustus.get_args(pipeline_args, genome).augustus_tm_gtf
                            for genome in tgt_genomes}
        elif mode == 'augTMR':
            tgt_genomes = pipeline_args.target_genomes
            gtf_in_files = {genome: Augustus.get_args(pipeline_args, genome).augustus_tmr_gtf
                            for genome in tgt_genomes}
        elif mode == 'transMap':
            tgt_genomes = pipeline_args.target_genomes
            gtf_in_files = {genome: FilterTransMap.get_args(pipeline_args, genome).filtered_tm_gtf
                            for genome in tgt_genomes}
        else:
            raise UserException('Invalid mode was passed to Hgm module: {}.'.format(mode))
        args = argparse.Namespace()
        args.genomes = tgt_genomes
        args.ref_genome = pipeline_args.ref_genome
        args.hal = pipeline_args.hal
        args.in_gtf = gtf_in_files
        args.gtf_out_dir = base_dir
        args.gtf_out_files = {genome: os.path.join(base_dir, genome + '.gtf') for genome in tgt_genomes}
        args.hints_db = pipeline_args.augustus_hints_db
        args.annotation_gp = ReferenceFiles.get_args(pipeline_args).annotation_gp
        # calculate the number of cores a hgm run should use
        # this is sort of a hack, but the reality is that halLiftover uses a fraction of a CPU most of the time
        max_cpu = min(pipeline_args.max_cores, multiprocessing.cpu_count())
        args.num_cpu = int(tools.mathOps.format_ratio(max_cpu, len(pipeline_args.modes)))
        return args

    def validate(self, pipeline_args):
        if not tools.misc.is_exec('homGeneMapping'):
            raise ToolMissingException('auxiliary program homGeneMapping from the Augustus package not in global path.')
        if not tools.misc.is_exec('halLiftover'):
            raise ToolMissingException('halLiftover from the halTools package not in global path.')
        if not tools.misc.is_exec('bedtools'):
            raise ToolMissingException('bedtools is required for the homGeneMapping module.')
        if pipeline_args.augustus_hints_db is None:
            raise InvalidInputException('Cannot run homGeneMapping module without a hints database.')
        if tools.hintsDatabaseInterface.hints_db_has_rnaseq(pipeline_args.augustus_hints_db) is False:
            raise UserException('homGeneMapping should not be ran on a hints database without RNA-seq.')

    def requires(self):
        pipeline_args = self.get_pipeline_args()
        self.validate(pipeline_args)
        for mode in pipeline_args.modes:
            yield self.clone(HgmDriverTask, mode=mode)


class HgmDriverTask(PipelineTask):
    """
    Task for running each individual instance of the Hgm pipeline. Dumps the results into a sqlite database
    Also produces a GTF file that is parsed into this database, but this file is not explicitly tracked by Luigi.
    """
    mode = luigi.Parameter()

    def output(self):
        pipeline_args = self.get_pipeline_args()
        hgm_args = Hgm.get_args(pipeline_args, self.mode)
        for genome in hgm_args.genomes:
            db = pipeline_args.dbs[genome]
            tools.fileOps.ensure_file_dir(db)
            conn_str = 'sqlite:///{}'.format(db)
            tablename = tools.sqlInterface.tables['hgm'][self.mode].__tablename__
            yield luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                      target_table=tablename,
                                                      update_id='_'.join([tablename, str(hash(pipeline_args))]))

    def requires(self):
        if self.no_evaluate_dependency is True:
            return
        if self.mode == 'augCGP':
            yield self.clone(AugustusCgp)
        elif self.mode == 'augTM' or self.mode == 'augTMR':
            yield self.clone(Augustus)
        elif self.mode == 'transMap':
            yield self.clone(FilterTransMap)
        else:
            raise UserException('Invalid mode passed to HgmDriverTask: {}.'.format(self.mode))

    def run(self):
        logger.info('Launching homGeneMapping for {}.'.format(self.mode))
        pipeline_args = self.get_pipeline_args()
        hgm_args = Hgm.get_args(pipeline_args, self.mode)
        hgm(hgm_args)
        # convert the output to a dataframe and write to the genome database
        databases = self.__class__.get_databases(pipeline_args)
        tablename = tools.sqlInterface.tables['hgm'][self.mode].__tablename__
        for genome, sqla_target in itertools.izip(*[hgm_args.genomes, self.output()]):
            # stores a dict with key=transcript_id and value=intron_support_count (e.g. "4,3,4,5")
            df = parse_hgm_gtf(hgm_args.gtf_out_files[genome])
            with tools.sqlite.ExclusiveSqlConnection(databases[genome]) as engine:
                df.to_sql(tablename, engine, if_exists='replace')
            sqla_target.touch()
            logger.info('Loaded table: {}.{}'.format(genome, tablename))


class AlignTranscripts(PipelineWrapperTask):
    """
    Aligns the transcripts from transMap/AugustusTMR/AugustusCGP to the parent transcript(s).
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'transcript_alignment')
        args = argparse.Namespace()
        args.ref_genome = pipeline_args.ref_genome
        args.genome = genome
        args.ref_genome_fasta = GenomeFiles.get_args(pipeline_args, pipeline_args.ref_genome).fasta
        args.genome_fasta = GenomeFiles.get_args(pipeline_args, genome).fasta
        args.annotation_gp = ReferenceFiles.get_args(pipeline_args).annotation_gp
        args.ref_db_path = PipelineTask.get_database(pipeline_args, pipeline_args.ref_genome)
        # the alignment_modes members hold the input genePreds and the mRNA/CDS alignment output paths
        args.transcript_modes = {'transMap': {'gp': FilterTransMap.get_args(pipeline_args, genome).filtered_tm_gp,
                                              'mRNA': os.path.join(base_dir, genome + '.transMap.mRNA.psl'),
                                              'CDS': os.path.join(base_dir, genome + '.transMap.CDS.psl')}}
        if pipeline_args.augustus is True:
            args.transcript_modes['augTM'] = {'gp':  Augustus.get_args(pipeline_args, genome).augustus_tm_gp,
                                              'mRNA': os.path.join(base_dir, genome + '.augTM.mRNA.psl'),
                                              'CDS': os.path.join(base_dir, genome + '.augTM.CDS.psl')}
        if pipeline_args.augustus_tmr is True:
            args.transcript_modes['augTMR'] = {'gp': Augustus.get_args(pipeline_args, genome).augustus_tmr_gp,
                                               'mRNA': os.path.join(base_dir, genome + '.augTMR.mRNA.psl'),
                                               'CDS': os.path.join(base_dir, genome + '.augTMR.CDS.psl')}
        if pipeline_args.augustus_cgp is True:
            args.transcript_modes['augCGP'] = {'gp': AugustusCgp.get_args(pipeline_args).augustus_cgp_gp[genome],
                                               'CDS': os.path.join(base_dir, genome + '.augCGP.CDS.psl')}
        return args

    def validate(self):
        if not tools.misc.is_exec('blat'):
            raise ToolMissingException('BLAT alignment tool not in global path.')

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            yield self.clone(AlignTranscriptDriverTask, genome=target_genome)


class AlignTranscriptDriverTask(ToilTask):
    """
    Task for per-genome launching of a toil pipeline for aligning all transcripts found back to the reference in
    transcript space using BLAT.

    Each task returns a PSL of all alignments that will be analyzed next by EvaluateTranscripts. For CGP transcripts,
    there may be more than one alignment.
    """
    genome = luigi.Parameter()

    def output(self):
        alignment_args = self.get_module_args(AlignTranscripts, genome=self.genome)
        for mode, paths in alignment_args.transcript_modes.iteritems():
            for aln_type in ['CDS', 'mRNA']:
                if aln_type in paths:
                    yield luigi.LocalTarget(paths[aln_type])

    def requires(self):
        if self.no_evaluate_dependency is True:
            return
        alignment_args = self.get_module_args(AlignTranscripts, genome=self.genome)
        if 'augTM' in alignment_args.transcript_modes:
            yield self.clone(Augustus)
        if 'augCGP' in alignment_args.transcript_modes:
            yield self.clone(AugustusCgp)
        yield self.clone(FilterTransMap)
        yield self.clone(ReferenceFiles)

    def run(self):
        logger.info('Launching Align Transcript toil pipeline for {} using {}.'.format(self.genome, self.batchSystem))
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'transcript_alignment', self.genome)
        toil_options = self.prepare_toil_options(toil_work_dir)
        alignment_args = self.get_module_args(AlignTranscripts, genome=self.genome)
        align_transcripts(alignment_args, toil_options)
        logger.info('Align Transcript toil pipeline for {} completed.'.format(self.genome))


class EvaluateTranscripts(PipelineWrapperTask):
    """
    Evaluates all transcripts for important features. See the classify.py module for details on how this works.

    Each task will generate a genome-specific sqlite database. See the classify.py docstring for details.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        args = argparse.Namespace()
        args.db_path = pipeline_args.dbs[genome]
        args.ref_db_path = PipelineTask.get_database(pipeline_args, pipeline_args.ref_genome)
        args.annotation_gp = ReferenceFiles.get_args(pipeline_args).annotation_gp
        args.fasta = GenomeFiles.get_args(pipeline_args, genome).fasta
        args.genome = genome
        args.ref_genome = pipeline_args.ref_genome
        # pass along all of the paths from alignment
        args.transcript_modes = AlignTranscripts.get_args(pipeline_args, genome).transcript_modes 
        return args

    def validate(self):
        pass

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        PipelineTask.get_database(pipeline_args, pipeline_args.ref_genome)
        for target_genome in pipeline_args.target_genomes:
            yield self.clone(EvaluateDriverTask, genome=target_genome)


class EvaluateDriverTask(PipelineTask):
    """
    Task for per-genome launching of a toil pipeline for aligning transcripts to their parent.
    """
    genome = luigi.Parameter()

    def build_table_names(self, eval_args):
        """construct table names based on input arguments"""
        tables = []
        for aln_mode in ['mRNA', 'CDS']:
            for tx_mode in eval_args.transcript_modes.iterkeys():
                if tx_mode == 'augCGP' and aln_mode == 'mRNA':
                    continue
                names = [x.__tablename__ for x in tools.sqlInterface.tables[aln_mode][tx_mode].values()]
                tables.extend(names)
        return tables

    def pair_table_output(self, eval_args):
        """return dict of {table_name: SQLAlchemyTarget} for final writing"""
        return dict(zip(*[self.build_table_names(eval_args), self.output()]))

    def write_to_sql(self, results, eval_args):
        """Load the results into the SQLite database"""
        with tools.sqlite.ExclusiveSqlConnection(eval_args.db_path) as engine:
            for table, target in self.pair_table_output(eval_args).iteritems():
                if table not in results:
                    continue
                df = results[table]
                df.to_sql(table, engine, if_exists='replace')
                target.touch()
                logger.info('Loaded table: {}.{}'.format(self.genome, table))

    def output(self):
        pipeline_args = self.get_pipeline_args()
        eval_args = self.get_module_args(EvaluateTranscripts, genome=self.genome)
        tools.fileOps.ensure_file_dir(eval_args.db_path)
        conn_str = 'sqlite:///{}'.format(eval_args.db_path)
        for table in self.build_table_names(eval_args):
            yield luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                      target_table=table,
                                                      update_id='_'.join([table, str(hash(pipeline_args))]))

    def requires(self):
        if self.no_evaluate_dependency is True:
            return
        return self.clone(AlignTranscripts), self.clone(ReferenceFiles)

    def run(self):
        logger.info('Evaluating transcript alignments for {}.'.format(self.genome))
        eval_args = self.get_module_args(EvaluateTranscripts, genome=self.genome)
        results = classify(eval_args)
        # results should be a dictionary of {table: dataframe}
        self.write_to_sql(results, eval_args)


class Consensus(PipelineWrapperTask):
    """
    Construct the consensus gene sets making use of the classification databases.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.out_dir, 'consensus_gene_set')
        # grab the genePred of every mode
        args = argparse.Namespace()
        args.cgp_num_exons = pipeline_args.cgp_num_exons
        args.cgp_splice_support = pipeline_args.cgp_splice_support
        gp_list = [TransMap.get_args(pipeline_args, genome).tm_gp]
        if pipeline_args.augustus is True:
            gp_list.append(Augustus.get_args(pipeline_args, genome).augustus_tm_gp)
        if pipeline_args.augustus_tmr is True:
            gp_list.append(Augustus.get_args(pipeline_args, genome).augustus_tmr_gp)
        if pipeline_args.augustus_cgp is True:
            gp_list.append(AugustusCgp.get_args(pipeline_args).augustus_cgp_gp[genome])
        args.gp_list = gp_list
        args.transcript_modes = AlignTranscripts.get_args(pipeline_args, genome).transcript_modes.keys()
        args.augustus_cgp = pipeline_args.augustus_cgp
        args.db_path = pipeline_args.dbs[genome]
        args.ref_db_path = PipelineTask.get_database(pipeline_args, pipeline_args.ref_genome)
        args.hints_db_has_rnaseq = pipeline_args.hints_db_has_rnaseq
        args.annotation_gp = ReferenceFiles.get_args(pipeline_args).annotation_gp
        args.consensus_gp = os.path.join(base_dir, genome + '.gp')
        args.consensus_gp_info = os.path.join(base_dir, genome + '.gp_info')
        args.consensus_gff3 = os.path.join(base_dir, genome + '.gff3')
        args.metrics_json = os.path.join(PipelineTask.get_metrics_dir(pipeline_args, genome), 'consensus.json')
        return args

    def validate(self):
        pass

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            yield self.clone(ConsensusDriverTask, genome=target_genome)


class ConsensusDriverTask(PipelineTask):
    """
    Driver task for performing consensus finding.
    """
    genome = luigi.Parameter()

    def output(self):
        consensus_args = self.get_module_args(Consensus, genome=self.genome)
        return luigi.LocalTarget(consensus_args.consensus_gp), luigi.LocalTarget(consensus_args.metrics_json)

    def requires(self):
        pipeline_args = self.get_pipeline_args()
        if self.no_evaluate_dependency is True:
            return
        yield self.clone(EvaluateTransMap)
        yield self.clone(EvaluateTranscripts)
        if pipeline_args.augustus_tmr is True:
            yield self.clone(Hgm)

    def run(self):
        consensus_args = self.get_module_args(Consensus, genome=self.genome)
        logger.info('Generating consensus gene set for {}.'.format(self.genome))
        consensus_gp, metrics_json = self.output()
        metrics_dict = generate_consensus(consensus_args, self.genome)
        PipelineTask.write_metrics(metrics_dict, metrics_json)


class Plots(PipelineTask):
    """
    Produce final analysis plots
    """
    @staticmethod
    def get_args(pipeline_args):
        base_dir = os.path.join(pipeline_args.out_dir, 'plots')
        ordered_genomes = tools.hal.build_genome_order(pipeline_args.hal, pipeline_args.ref_genome,
                                                       genome_subset=pipeline_args.target_genomes)
        args = argparse.Namespace()
        args.ordered_genomes = ordered_genomes
        # plots derived from transMap results
        args.tm_coverage = luigi.LocalTarget(os.path.join(base_dir, 'transmap_coverage.pdf'))
        args.tm_identity = luigi.LocalTarget(os.path.join(base_dir, 'transmap_identity.pdf'))
        # plots derived from transMap filtering
        args.paralogy = luigi.LocalTarget(os.path.join(base_dir, 'paralogy.pdf'))
        args.transmap_filtering = luigi.LocalTarget(os.path.join(base_dir, 'transmap_filtering.pdf'))
        # plots derived from transcript alignment / consensus finding
        args.coverage = luigi.LocalTarget(os.path.join(base_dir, 'coverage.pdf'))
        args.identity = luigi.LocalTarget(os.path.join(base_dir, 'identity.pdf'))
        args.consensus_score = luigi.LocalTarget(os.path.join(base_dir, 'consensus_score.pdf'))
        args.completeness = luigi.LocalTarget(os.path.join(base_dir, 'completeness.pdf'))
        args.categories = luigi.LocalTarget(os.path.join(base_dir, 'transcript_categories.pdf'))
        args.gene_failure = luigi.LocalTarget(os.path.join(base_dir, 'gene_failure.pdf'))
        args.transcript_failure = luigi.LocalTarget(os.path.join(base_dir, 'transcript_failure.pdf'))
        args.consensus_splice_support = luigi.LocalTarget(os.path.join(base_dir, 'consensus_splice_support.pdf'))
        if 'augTM' in pipeline_args.modes or 'augTMR' in pipeline_args.modes or 'augCGP' in pipeline_args.modes:
            args.tx_modes = luigi.LocalTarget(os.path.join(base_dir, 'transcript_modes.pdf'))
        # plots that depend on execution mode
        if 'augCGP' in pipeline_args.modes:
            args.novel = luigi.LocalTarget(os.path.join(base_dir, 'cgp_novel.pdf'))
        if pipeline_args.resolve_split_genes is True:
            args.split_genes = luigi.LocalTarget(os.path.join(base_dir, 'split_genes.pdf'))
        # input data
        args.metrics_jsons = OrderedDict([[genome, Consensus.get_args(pipeline_args, genome).metrics_json]
                                          for genome in ordered_genomes])
        args.tm_jsons = OrderedDict([[genome, FilterTransMap.get_args(pipeline_args, genome).metrics_json]
                                     for genome in ordered_genomes])
        args.annotation_db = PipelineTask.get_database(pipeline_args, pipeline_args.ref_genome)
        args.dbs = OrderedDict([[genome, PipelineTask.get_database(pipeline_args, genome)]
                                for genome in ordered_genomes])
        return args

    def output(self):
        pipeline_args = self.get_pipeline_args()
        args = Plots.get_args(pipeline_args)
        return [p for p in args.__dict__.itervalues() if isinstance(p, luigi.LocalTarget)]

    def requires(self):
        yield self.clone(Consensus)

    def run(self):
        pipeline_args = self.get_pipeline_args()
        logger.info('Generating plots.')
        generate_plots(Plots.get_args(pipeline_args))


###
# Entry point without using luigi
###


def parse_args():
    """
    Provides the ability to run directly from this script, bypassing the luigi wrapper. If you go this route, you
    cannot control the number of concurrent toil pipelines.
    :return: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--hal', required=True)
    parser.add_argument('--ref-genome', required=True)
    parser.add_argument('--annotation', required=True)
    parser.add_argument('--out-dir', default='./cat_output')
    parser.add_argument('--work-dir', default=os.path.join(tempfile.gettempdir(), __name__))
    parser.add_argument('--target-genomes', nargs='+', default=None)
    # parallelism
    parser.add_argument('--workers', default=10, type=int)
    # Debugging option - use this to bypass the dependency graph for specific submodules
    parser.add_argument('--no-evaluate-dependency', action='store_true')
    # augustus TM(R) options
    parser.add_argument('--augustus', action='store_true')
    parser.add_argument('--augustus-species', default='human')
    parser.add_argument('--augustus-hints-db', default=None)
    parser.add_argument('--tm-cfg', default='augustus_cfgs/extrinsic.ETM1.cfg')
    parser.add_argument('--tmr-cfg', default='augustus_cfgs/extrinsic.ETM2.cfg')
    # augustus CGP options
    parser.add_argument('--augustus-cgp', action='store_true')
    parser.add_argument('--augustus-cgp-cfg-template', default='augustus_cfgs/cgp_extrinsic_template.cfg')
    parser.add_argument('--cgp-param', default='augustus_cfgs/log_reg_parameters_default.cfg')
    parser.add_argument('--maf-chunksize', default=2500000, type=int)
    parser.add_argument('--maf-overlap', default=500000, type=int)
    # consensus options
    parser.add_argument('--resolve-split-genes', action='store_true')
    parser.add_argument('--cgp-splice-support', default=0.8, type=float)
    parser.add_argument('--cgp-num-exons', default=3, type=int)
    # toil options
    parser.add_argument('--batchSystem', default='singleMachine')
    parser.add_argument('--maxCores', default=16, type=int)
    parser.add_argument('--logLevel', default='WARNING')
    parser.add_argument('--cleanWorkDir', default='onSuccess')
    parser.add_argument('--parasolCommand', default=None)
    args = parser.parse_args()
    args.target_genomes = tuple(args.target_genomes) if args.target_genomes is not None else None
    return args


if __name__ == '__main__':
    args = parse_args()
    workers = args.workers
    del args.workers  # hack because workers is actually a argument to luigi, not RunCat
    luigi.build([RunCat(**vars(args))], logging_conf_file='logging.cfg', workers=workers)
