#!/usr/bin/env python2.7
"""
Comparative Annotation Toolkit.
"""
import logging
import argparse
import itertools
import os
import luigi
import sqlalchemy
import sqlite3 as sql
import tempfile
import tools.fileOps
import tools.procOps
import tools.toilInterface
import tools.nameConversions
import tools.hal
import tools.transcripts
import tools.sqlInterface
import tools.psl
import tools.bio
import tools.gff3
import tools.misc
from tools.luigiAddons import multiple_inherits, multiple_requires, AbstractAtomicFileTask
from luigi.util import inherits, requires
import luigi.contrib.sqla
from chaining import chaining
from augustus import augustus
from augustus_cgp import augustus_cgp
from hgm import hgm, parse_hgm_gtf
from align_transcripts import align_transcripts
from classify import classify


logger = logging.getLogger(__name__)


class UserException(Exception):
    pass


class PipelineParameterMixin(object):
    """
    Mixin class for adding extracting the base arguments to the pipeline. Expects that the mixed in class is a 
    derivative of RunCat.
    
    This can't be a classmethod of RunCat due to how Luigi resolves parameters.
    """
    def get_pipeline_args(self):
        """returns a args of all of the arguments to the pipeline. Resolves the target genomes variable"""
        args = argparse.Namespace()
        args.hal = os.path.abspath(self.hal)
        args.ref_genome = self.ref_genome
        args.annotation = os.path.abspath(self.annotation)
        args.out_dir = os.path.abspath(self.out_dir)
        args.work_dir = os.path.abspath(self.work_dir)
        args.augustus = self.augustus
        args.augustus_cgp = self.augustus_cgp
        args.augustus_species = self.augustus_species
        args.tm_to_hints_script = self.tm_to_hints_script
        if self.augustus_hints_db is not None:
            args.augustus_hints_db = os.path.abspath(self.augustus_hints_db)
        else:
            args.augustus_hints_db = None
        args.tm_cfg = os.path.abspath(self.tm_cfg)
        args.tmr_cfg = os.path.abspath(self.tmr_cfg)
        args.tm_to_hints_script = os.path.abspath(self.tm_to_hints_script)
        args.augustus_cgp = self.augustus_cgp
        args.maf_chunksize = self.maf_chunksize
        args.maf_overlap = self.maf_overlap
        if self.augustus_cgp_cfg is not None:
            args.augustus_cgp_cfg = os.path.abspath(self.augustus_cgp_cfg)
        else:
            args.augustus_cgp_cfg = None
        if self.augustus_cgp_param is not None:
            args.augustus_cgp_param = os.path.abspath(self.augustus_cgp_param)
        else:
            args.augustus_cgp_param = None
        if self.target_genomes is None:
            target_genomes = tools.hal.extract_genomes(self.hal)
            target_genomes = tuple(set(target_genomes) - {self.ref_genome})
        else:
            target_genomes = tuple([x for x in self.target_genomes])
        args.target_genomes = target_genomes
        return args


class RunCat(luigi.WrapperTask, tools.toilInterface.ToilOptionsMixin):
    """
    Task that executes the entire pipeline.
    """
    hal = luigi.Parameter()
    ref_genome = luigi.Parameter()
    annotation = luigi.Parameter()
    out_dir = luigi.Parameter(default='./cat_output')
    work_dir = luigi.Parameter(default=tempfile.gettempdir())
    target_genomes = luigi.TupleParameter(default=None)
    # AugustusTM(R) parameters
    augustus = luigi.BoolParameter(default=False)
    augustus_species = luigi.Parameter(default='human')
    augustus_hints_db = luigi.Parameter(default=None)
    tm_cfg = luigi.Parameter(default='augustus_cfgs/extrinsic.ETM1.cfg', significant=False)
    tmr_cfg = luigi.Parameter(default='augustus_cfgs/extrinsic.ETM2.cfg', significant=False)
    # AugustusCGP parameters
    augustus_cgp = luigi.BoolParameter(default=False)
    augustus_cgp_cfg = luigi.Parameter(default=None, significant=False)
    augustus_cgp_param = luigi.Parameter(default='augustus_cfgs/log_reg_parameters_default.cfg', significant=False)
    maf_chunksize = luigi.IntParameter(default=2500000, significant=False)
    maf_overlap = luigi.IntParameter(default=500000, significant=False)
    tm_to_hints_script = luigi.Parameter(default='tools/transMap2hints.pl', significant=False)

    def requires(self):
        yield self.clone(PrepareFiles)
        yield self.clone(Chaining)
        yield self.clone(TransMap)
        if self.augustus is True:
            yield self.clone(Augustus)
        if self.augustus_cgp is True:
            yield self.clone(AugustusCgp)
        yield self.clone(Hgm)
        yield self.clone(AlignTranscripts)
        yield self.clone(EvaluateTranscripts)
        #yield self.clone(Consensus)
        #yield self.clone(Plots)


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
class GenomeFiles(luigi.WrapperTask, PipelineParameterMixin):
    """
    WrapperTask for producing all genome files.

    GenomeFiles -> GenomeFasta -> GenomeTwoBit -> GenomeFlatFasta -> GenomeFastaIndex
                -> GenomeSizes

    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'genome_files')
        args = {'fasta': os.path.join(base_dir, genome + '.fa'),
                'two_bit': os.path.join(base_dir, genome + '.2bit'),
                'fasta_index': os.path.join(base_dir, genome + '.fa.fai'),
                'sizes': os.path.join(base_dir, genome + '.chrom.sizes'),
                'flat_fasta': os.path.join(base_dir, genome + '.fa.flat'),
                'genome': genome}
        return args

    def requires(self):
        pipeline_args = self.get_pipeline_args()
        for genome in itertools.chain(pipeline_args.target_genomes, [self.ref_genome]):
            args = self.get_args(pipeline_args, genome)
            yield self.clone(GenomeFasta, **args)
            yield self.clone(GenomeTwoBit, **args)
            yield self.clone(GenomeSizes, **args)
            yield self.clone(GenomeFlatFasta, **args)


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


@inherits(PrepareFiles)
class ReferenceFiles(luigi.WrapperTask, PipelineParameterMixin):
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
        args = {'annotation_gp': os.path.join(base_dir, annotation + '.gp'),
                'annotation_db': os.path.join(base_dir, annotation + '.db'),
                'transcript_fasta': os.path.join(base_dir, annotation + '.fa'),
                'transcript_flat_fasta': os.path.join(base_dir, annotation + '.fa.flat'),
                'transcript_bed': os.path.join(base_dir, annotation + '.bed'),
                'ref_psl': os.path.join(base_dir, annotation + '.psl')}
        args.update(GenomeFiles.get_args(pipeline_args, pipeline_args.ref_genome))
        return args

    def requires(self):
        pipeline_args = self.get_pipeline_args()
        args = self.get_args(pipeline_args)
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
        logger.info('Converting annotation gff3 to genePred.')
        cmd = ['gff3ToGenePred', '-rnaNameAttr=transcript_id', '-geneNameAttr=gene_id', '-honorStartStopCodons',
               self.annotation, '/dev/stdout']
        self.run_cmd(cmd)


@inherits(ReferenceFiles)
class Gff3ToAttrs(luigi.Task):
    """
    Uses the gff3 parser to extract the attributes table, converting the table into a sqlite databse for import
    and use by downstream tools.
    """
    annotation_db = luigi.Parameter()

    def output(self):
        tools.fileOps.ensure_file_dir(self.annotation_db)
        conn_str = 'sqlite:///{}'.format(self.annotation_db)
        attrs_table = luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                          target_table=self.ref_genome,
                                                          update_id=self.task_id)
        return attrs_table

    def run(self):
        logger.info('Extracting gff3 attributes to sqlite database.')
        results = tools.gff3.extract_attrs(self.annotation)
        engine = sqlalchemy.create_engine('sqlite:///{}'.format(self.annotation_db))
        results.to_sql(self.ref_genome, engine, if_exists='replace')
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
    Produces a fasta for each transcript. Requires bedtools.

    TODO: bedtools 2.26 broke this. -name now does not just provide the sequence name. I have to use my transcript
    library.
    """
    transcript_fasta = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(self.transcript_fasta)

    def run(self):
        logger.info('Extracting reference annotation fasta.')
        #cmd = ['bedtools', 'getfasta', '-fi', self.fasta, '-bed', self.transcript_bed, '-fo', '/dev/stdout',
        #       '-name', '-split', '-s']
        #self.run_cmd(cmd)
        seq_dict = tools.bio.get_sequence_dict(self.fasta, upper=False)
        seqs = {name: tx.get_mrna(seq_dict) for name, tx in tools.transcripts.transcript_iterator(self.transcript_bed)}
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


@inherits(RunCat)
class Chaining(luigi.WrapperTask, PipelineParameterMixin):
    """
    WrapperTask for producing chain files for each combination of ref_genome-target_genome
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'chaining')
        args = {'hal': pipeline_args.hal,
                'ref_genome': pipeline_args.ref_genome, 'genome': genome,
                'chain_file': os.path.join(base_dir, '{}-{}.chain'.format(pipeline_args.ref_genome, genome))}
        ref_files = GenomeFiles.get_args(pipeline_args, pipeline_args.ref_genome)
        tgt_files = GenomeFiles.get_args(pipeline_args, genome)
        args['query_two_bit'] = ref_files['two_bit']
        args['target_two_bit'] = tgt_files['two_bit']
        args['query_sizes'] = ref_files['sizes']
        return args

    def validate(self):
        # TODO: make sure that all input args exist
        pass

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            args = self.get_args(pipeline_args, target_genome)
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
        logger.info('Launching Pairwise Chaining toil pipeline for {}.'.format(self.genome))
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'chaining', self.genome)
        toil_options = self.prepare_toil_options(toil_work_dir)
        chaining(self.chain_args, toil_options)
        logger.info('Pairwise Chaining for {} is complete.'.format(self.genome))


@inherits(RunCat)
class TransMap(luigi.WrapperTask, PipelineParameterMixin):
    """
    Runs transMap.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        base_dir = os.path.join(pipeline_args.work_dir, 'transMap')
        ref_files = ReferenceFiles.get_args(pipeline_args)
        tgt_genome_files = GenomeFiles.get_args(pipeline_args, genome)
        chain_args = Chaining.get_args(pipeline_args, genome)
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
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            args = self.get_args(pipeline_args, target_genome)
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
        logger.info('Running transMap for {}.'.format(self.genome))
        tools.fileOps.ensure_file_dir(self.output().path)
        psl_cmd = ['pslMap', '-chainMapFile', self.tm_args['ref_psl'],
                   self.tm_args['chain_file'], '/dev/stdout']
        post_chain_cmd = ['postTransMapChain', '/dev/stdin', '/dev/stdout']
        sort_cmd = ['sort', '-k', '14,14', '-k', '16,16n']
        recalc_cmd = ['pslRecalcMatch', '/dev/stdin', self.tm_args['two_bit'], self.tm_args['transcript_fasta'],
                      'stdout']
        cmd_list = [psl_cmd, post_chain_cmd, sort_cmd, recalc_cmd]
        # hacky way to make unique - capture output to a file, then process
        logger.info('Writing transMap results to file for {}.'.format(self.genome))
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
        logger.info('Converting transMap PSL to genePred for {}.'.format(self.genome))
        cmd = ['transMapPslToGenePred', '-nonCodingGapFillMax=80', '-codingGapFillMax=50',
               self.tm_args['annotation_gp'], self.tm_args['tm_psl'], '/dev/stdout']
        self.run_cmd(cmd)


@inherits(RunCat)
class Augustus(luigi.WrapperTask, PipelineParameterMixin):
    """
    Runs AugustusTM(R) on the coding output from transMap.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        tm_args = TransMap.get_args(pipeline_args, genome)
        tgt_genome_files = GenomeFiles.get_args(pipeline_args, genome)
        annotation_files = ReferenceFiles.get_args(pipeline_args)
        base_dir = os.path.join(pipeline_args.work_dir, 'augustus')
        args = {'ref_genome': pipeline_args.ref_genome, 'genome': genome,
                'genome_fasta': tgt_genome_files['fasta'],
                'ref_psl': tm_args['ref_psl'],
                'annotation_gp': annotation_files['annotation_gp'],
                'annotation_db': annotation_files['annotation_db'],
                'tm_gp': tm_args['tm_gp'],
                'tm_psl': tm_args['tm_psl'],
                'augustus_tm_gp': os.path.join(base_dir, genome + '.TM.gp'),
                'augustus_tm_gtf': os.path.join(base_dir, genome + '.TM.gtf'),
                'augustus_hints_db': pipeline_args.augustus_hints_db if pipeline_args.augustus_hints_db is not None else None,
                'tm_cfg': pipeline_args.tm_cfg,
                'tmr_cfg': pipeline_args.tmr_cfg,
                'tm_to_hints_script': pipeline_args.tm_to_hints_script,
                'augustus_species': pipeline_args.augustus_species}
        if pipeline_args.augustus_hints_db is not None:
            args['augustus_tmr_gp'] = os.path.join(base_dir, genome + '.TMR.gp')
            args['augustus_tmr_gtf'] = os.path.join(base_dir, genome + '.TMR.gtf')
        return args

    def validate(self):
        # TODO: make sure that all input args exist
        pass

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            args = self.get_args(pipeline_args, target_genome)
            yield self.clone(AugustusDriverTask, augustus_args=args, genome=target_genome)


@inherits(Augustus)
class AugustusDriverTask(tools.toilInterface.ToilTask):
    """
    Task for per-genome launching of a toil pipeline for running Augustus.
    """
    augustus_args = luigi.DictParameter()  # dict to pass directly to TMR toil module
    genome = luigi.Parameter()

    def output(self):
        yield luigi.LocalTarget(self.augustus_args['augustus_tm_gp'])
        yield luigi.LocalTarget(self.augustus_args['augustus_tm_gtf'])
        if 'augustus_tmr_gp' in self.augustus_args:
            yield luigi.LocalTarget(self.augustus_args['augustus_tmr_gp'])
            yield luigi.LocalTarget(self.augustus_args['augustus_tmr_gtf'])

    def requires(self):
        return self.clone(TransMap)

    def extract_coding_genes(self):
        """extracts only coding genes from the input genePred, returning a path to a tmp file"""
        coding_gp = tools.fileOps.get_tmp_file()
        attrs = tools.sqlInterface.read_attrs(self.augustus_args['annotation_db'], self.augustus_args['ref_genome'])
        names = set(attrs[attrs.tx_biotype == 'protein_coding'].index)
        with open(coding_gp, 'w') as outf:
            for name, tx in tools.transcripts.gene_pred_iterator(self.augustus_args['tm_gp']):
                if tools.nameConversions.strip_alignment_numbers(name) in names:
                    tools.fileOps.print_row(outf, tx.get_gene_pred())
        if os.path.getsize(coding_gp) == 0:
            raise RuntimeError('Unable to extract coding transcripts from the transMap genePred.')
        return coding_gp

    def run(self):
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'augustus', self.genome)
        logger.info('Launching Augustus toil pipeline on {}.'.format(self.genome))
        toil_options = self.prepare_toil_options(toil_work_dir)
        coding_gp = self.extract_coding_genes()
        augustus(self.augustus_args, coding_gp, toil_options)
        logger.info('Augustus toil pipeline for {} completed.'.format(self.genome))
        os.remove(coding_gp)
        for out_gp, out_gtf in tools.misc.pairwise(self.output()):
            tools.misc.convert_gtf_gp(out_gp, out_gtf)


@inherits(RunCat)
class AugustusCgp(tools.toilInterface.ToilTask, PipelineParameterMixin):
    """
    Task for launching the AugustusCGP toil pipeline
    """
    @staticmethod
    def get_args(pipeline_args):
        # add reference to the target genomes
        tgt_genomes = list(pipeline_args.target_genomes) + [pipeline_args.ref_genome]
        fasta_files = {genome: GenomeFiles.get_args(pipeline_args, genome)['fasta'] for genome in tgt_genomes}
        base_dir = os.path.join(pipeline_args.work_dir, 'augustus_cgp')
        # output
        output_gp_files = {genome: os.path.abspath(os.path.join(base_dir, genome + '_augustus_cgp.gp'))
                           for genome in tgt_genomes}
        gtf_files = {genome: os.path.abspath(os.path.join(base_dir, genome + '_augustus_cgp.gtf'))
                     for genome in tgt_genomes}
        # transMap files used for assigning parental gene
        tm_gp_files = {genome: TransMap.get_args(pipeline_args, genome)['tm_gp']
                       for genome in pipeline_args.target_genomes}
        ref_files = ReferenceFiles.get_args(pipeline_args)
        # add the reference annotation as a pseudo-transMap to assign parents in reference
        tm_gp_files[pipeline_args.ref_genome] = ref_files['annotation_gp']
        cgp_cfg = os.path.abspath(pipeline_args.augustus_cgp_cfg) if pipeline_args.augustus_cgp_cfg is not None else None
        hints_db = os.path.abspath(pipeline_args.augustus_hints_db) if pipeline_args.augustus_hints_db is not None else None
        args = {'genomes': tgt_genomes,
                'fasta_files': fasta_files,
                'tm_gps': tm_gp_files,
                'hal': pipeline_args.hal,
                'ref_genome': pipeline_args.ref_genome,
                'augustus_cgp_gp': output_gp_files,
                'augustus_cgp_gtf': gtf_files,
                'species': pipeline_args.augustus_species,
                'chunksize': pipeline_args.maf_chunksize,
                'overlap': pipeline_args.maf_overlap,
                'cgp_cfg': cgp_cfg,
                'hints_db': hints_db,
                'cgp_param': os.path.abspath(pipeline_args.augustus_cgp_param),
                'annotation_db': ref_files['annotation_db']}
        # chromSizes: required for splitting the HAL alignment along the reference genome
        ref_file = GenomeFiles.get_args(pipeline_args, pipeline_args.ref_genome)
        args['query_sizes'] = ref_file['sizes']
        return args

    def output(self):
        pipeline_args = self.get_pipeline_args()
        cgp_args = self.get_args(pipeline_args)
        for cat in ['augustus_cgp_gp', 'augustus_cgp_gtf']:
            for path in cgp_args[cat].itervalues():
                yield luigi.LocalTarget(path)

    def validate(self):
        # make sure that all external tools are executable
        if not tools.misc.is_exec('joingenes'):
            raise UserException('auxiliary program joingenes from the Augustus package not in global path.')
        if not tools.misc.is_exec('augustus'):
            raise UserException('augustus not in global path.')
        if not tools.misc.is_exec('hal2maf'):
            raise UserException('hal2maf from the halTools package not in global path.')
        # TODO: make sure that all input args exist

    def requires(self):
        self.validate()
        yield self.clone(TransMap)

    def run(self):
        pipeline_args = self.get_pipeline_args()
        if pipeline_args.augustus_hints_db is None:
            # TODO: if no database is given (de novo setting),
            # create a new DB with the flattened genomes from the HAL alignment
            raise UserException('Cannot run AugustusCGP without a hints database.')
        logger.info('Beginning AugustusCGP toil pipeline.')
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'augustus_cgp')
        toil_options = self.prepare_toil_options(toil_work_dir)
        cgp_args = self.get_args(pipeline_args)
        augustus_cgp(cgp_args, toil_options)
        logger.info('AugustusCGP toil pipeline completed.')
        # convert each to genePred as well
        for genome in itertools.chain(pipeline_args.target_genomes, [pipeline_args.ref_genome]):
            gp_target = luigi.LocalTarget(cgp_args['augustus_cgp_gp'][genome])
            gtf_target = luigi.LocalTarget(cgp_args['augustus_cgp_gtf'][genome])
            tools.misc.convert_gtf_gp(gp_target, gtf_target)
        logger.info('Finished converting AugustusCGP output.')


@inherits(RunCat)
class Hgm(tools.toilInterface.ToilTask, PipelineParameterMixin):
    """
    Task for launching the HomGeneMapping toil pipeline
    """
    @staticmethod
    def get_args(pipeline_args):
        # add reference to the target genomes
        tgt_genomes = list(pipeline_args.target_genomes) + [pipeline_args.ref_genome]
        # input (needs to be changed to the concatenated gene sets from AugustusTM(R), AugustusCgp (and TransMap?))
        gtf_in_files = {genome: AugustusCgp.get_args(pipeline_args)['augustus_cgp_gtf'][genome]
                     for genome in tgt_genomes}
        base_dir = os.path.join(pipeline_args.work_dir, 'hgm')
        # output (only temporary for debugging)
        gtf_out_files = {genome: os.path.abspath(os.path.join(base_dir, genome + '_hgm.gtf'))
                     for genome in tgt_genomes}
        gp_out_files = {genome: os.path.abspath(os.path.join(base_dir, genome + '_hgm.gp'))
                     for genome in tgt_genomes}
        hints_db = os.path.abspath(pipeline_args.augustus_hints_db) if pipeline_args.augustus_hints_db is not None else None
        args = {'genomes': tgt_genomes,
                'hal': pipeline_args.hal,
                'hgm_gtf': gtf_out_files,
                'hgm_gp': gp_out_files,
                'in_gtf': gtf_in_files,
                'hints_db': hints_db}
        return args

    def output(self):
        pipeline_args = self.get_pipeline_args()
        hgm_args = self.get_args(pipeline_args)
        for cat in ['hgm_gp', 'hgm_gtf']:
            for path in hgm_args[cat].itervalues():
                yield luigi.LocalTarget(path)

    def validate(self):
        # make sure that all external tools are executable
        if not tools.misc.is_exec('homGeneMapping'):
            raise UserException('auxiliary program homGeneMapping from the Augustus package not in global path.')
        if not tools.misc.is_exec('halLiftover'):
            raise UserException('halLiftover from the halTools package not in global path.')

    def requires(self):
        self.validate()
        yield self.clone(AugustusCgp) # add AugustusTM(R)/TransMap if necessary

    def appendSJSupportToGp(self, gp, intron_support_counts):
        """
        appends a column with the intron support counts to a genePred file
        """
        out_lines = []
        with open(gp) as inf:
                for line in inf:
                    txid = line.split()[0]
                    counts = intron_support_counts.get(txid, "") # empty string for single exon genes
                    out_lines.append("\t".join([line.strip(),counts]))
        with open(gp, "w") as outf:
            outf.write("\n".join(out_lines))

    def run(self):
        pipeline_args = self.get_pipeline_args()
        logger.info('Beginning Hgm toil pipeline.')
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'hgm')
        toil_options = self.prepare_toil_options(toil_work_dir)
        hgm_args = self.get_args(pipeline_args)
        hgm(hgm_args, toil_options)
        logger.info('Hgm toil pipeline completed.')
        # this will be the 
        intron_support_counts = {genome: parse_hgm_gtf(path) for genome,path in hgm_args['hgm_gtf'].iteritems()}
        # convert input gtf to gp and add a column with the intron support counts to th gp (only temporary for debugging)
        for genome in itertools.chain(hgm_args['genomes']):
            gtf = luigi.LocalTarget(hgm_args['in_gtf'][genome])
            gp = luigi.LocalTarget(hgm_args['hgm_gp'][genome])
            tools.misc.convert_gtf_gp(gp, gtf)
            self.appendSJSupportToGp(hgm_args['hgm_gp'][genome], intron_support_counts[genome])
        logger.info('Finished converting Hgm output.')   


@inherits(RunCat)
class AlignTranscripts(luigi.WrapperTask, PipelineParameterMixin):
    """
    Aligns the transcripts from transMap/AugustusTMR/AugustusCGP to the parent transcript(s).
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        tgt_genome_files = GenomeFiles.get_args(pipeline_args, genome)
        ref_genome_files = GenomeFiles.get_args(pipeline_args, pipeline_args.ref_genome)
        tm_files = TransMap.get_args(pipeline_args, genome)
        annotation_files = ReferenceFiles.get_args(pipeline_args)
        base_dir = os.path.join(pipeline_args.work_dir, 'transcript_alignment')
        args = {'ref_genome': pipeline_args.ref_genome, 'genome': genome,
                'ref_genome_fasta': ref_genome_files['fasta'],
                'annotation_gp': annotation_files['annotation_gp'],
                'annotation_db': annotation_files['annotation_db'],
                'genome_fasta': tgt_genome_files['fasta'],
                'modes': {'transMap': {'gp': tm_files['tm_gp'],
                                       'MUSCLE': os.path.join(base_dir, genome + '.transMap.MUSCLE.fasta.gz'),
                                       'PRANK': os.path.join(base_dir, genome + '.transMap.PRANK.fasta.gz')}}
                }
        if pipeline_args.augustus is True:
            aug_args = Augustus.get_args(pipeline_args, genome)
            args['modes']['augTM'] = {'gp': aug_args['augustus_tm_gp'],
                                      'MUSCLE': os.path.join(base_dir, genome + '.augTM.MUSCLE.fasta.gz'),
                                      'PRANK': os.path.join(base_dir, genome + '.augTM.PRANK.fasta.gz')}
            if pipeline_args.augustus_hints_db is not None:
                args['modes']['augTMR'] = {'gp': aug_args['augustus_tmr_gp'],
                                           'MUSCLE': os.path.join(base_dir, genome + '.augTMR.MUSCLE.fasta.gz'),
                                           'PRANK': os.path.join(base_dir, genome + '.augTMR.PRANK.fasta.gz')}
        if pipeline_args.augustus_cgp is True:
            cgp_args = AugustusCgp.get_args(pipeline_args)
            args['modes']['augCGP'] = {'gp': cgp_args['augustus_cgp_gp'][genome],
                                       'PRANK': os.path.join(base_dir, genome + '.augCGP.PRANK.fasta.gz')}
        return args

    def validate(self):
        # TODO: make sure that all input args exist
        pass

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            args = self.get_args(pipeline_args, target_genome)
            yield self.clone(AlignTranscriptDriverTask, alignment_args=args, genome=target_genome)


@inherits(AlignTranscripts)
class AlignTranscriptDriverTask(tools.toilInterface.ToilTask):
    """
    Task for per-genome launching of a toil pipeline for aligning all transcripts found back to the reference in
    transcript space using BLAT.

    Each task returns a PSL of all alignments that will be analyzed next by EvaluateTranscripts. For CGP transcripts,
    there may be more than one alignment.
    """
    alignment_args = luigi.DictParameter()  # dict to pass directly to alignment toil module
    genome = luigi.Parameter()

    def output(self):
        for mode, paths in self.alignment_args['modes'].iteritems():
            for aln_type in ['PRANK', 'MUSCLE']:
                if aln_type in paths:
                    yield luigi.LocalTarget(paths[aln_type])

    def requires(self):
        if 'augTM' in self.alignment_args['modes']:
            yield self.clone(Augustus)
        if 'augCGP' in self.alignment_args['modes']:
            yield self.clone(AugustusCgp)
        yield self.clone(TransMap)

    def run(self):
        logger.info('Beginning Align Transcript toil pipeline for {}.'.format(self.genome))
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'transcript_alignment', self.genome)
        toil_options = self.prepare_toil_options(toil_work_dir)
        align_transcripts(self.alignment_args, toil_options)
        logger.info('Align Transcript toil pipeline for {} completed.'.format(self.genome))


@inherits(RunCat)
class EvaluateTranscripts(luigi.WrapperTask, PipelineParameterMixin):
    """
    Evaluates all transcripts for important features. See the classify.py module for details on how this works.

    Each task will generate a genome-specific sqlite database. See the classify.py docstring for details.
    """
    @staticmethod
    def get_args(pipeline_args, genome):
        tm_args = TransMap.get_args(pipeline_args, genome)
        tgt_genome_files = GenomeFiles.get_args(pipeline_args, genome)
        annotation_files = ReferenceFiles.get_args(pipeline_args)
        tx_alignment_args = AlignTranscripts.get_args(pipeline_args, genome)
        base_out_dir = os.path.join(pipeline_args.out_dir, 'databases')
        args = {'db': os.path.join(base_out_dir, '{}.db'.format(genome)),
                'tm_psl': tm_args['tm_psl'],
                'tm_gp': tm_args['tm_gp'],
                'annotation_gp': annotation_files['annotation_gp'],
                'annotation_db': annotation_files['annotation_db'],
                'genome_fasta': tgt_genome_files['fasta'],
                'genome': genome,
                'ref_genome': pipeline_args.ref_genome,
                'modes': tx_alignment_args['modes']}  # pass along all of the paths from evaluation
        return args

    def validate(self):
        # TODO: make sure that all input args exist
        pass

    def requires(self):
        self.validate()
        pipeline_args = self.get_pipeline_args()
        for target_genome in pipeline_args.target_genomes:
            args = self.get_args(pipeline_args, target_genome)
            yield self.clone(EvaluateDriverTask, eval_args=args, genome=target_genome)


@multiple_inherits(EvaluateTranscripts)
class EvaluateDriverTask(tools.toilInterface.ToilTask):
    """
    Task for per-genome launching of a toil pipeline for aligning transcripts to their parent.
    """
    eval_args = luigi.DictParameter()  # dict to pass directly to evaluate toil module
    genome = luigi.Parameter()

    def build_table_names(self):
        """construct table names based on input arguments"""
        tables = ['Alignment']
        for tx_mode, path_dict in self.eval_args['modes'].iteritems():
            for aln_mode in ['MUSCLE', 'PRANK']:
                if aln_mode in path_dict:
                    metrics_table = '_'.join([aln_mode, tx_mode, 'Metrics'])
                    eval_table = '_'.join([aln_mode, tx_mode, 'Evaluation'])
                    tables.extend([metrics_table, eval_table])
        return tables

    def pair_table_output(self):
        """return dict of {table_name: SQLAlchemyTarget} for final writing"""
        return dict(zip(*[self.build_table_names(), self.output()]))

    def write_to_sql(self, results):
        """Load the results into the SQLite database"""
        # TODO: sqlalchemy can't write hierarchical indices. See this:
        # http://stackoverflow.com/questions/39089513/pandas-sqlalchemy-engine-does-not-produce-hierarchical-index-but-legacy-mode-do
        #engine = sqlalchemy.create_engine('sqlite:///{}'.format(self.eval_args['db']))
        engine = sql.connect(self.eval_args['db'])
        for table, target in self.pair_table_output().iteritems():
            if table not in results:
                continue
            logger.info('Loading table: {}.{}'.format(self.genome, table))
            df = results[table]
            df.to_sql(table, engine, if_exists='replace')
            target.touch()
        engine.close()

    def output(self):
        tools.fileOps.ensure_file_dir(self.eval_args['db'])
        conn_str = 'sqlite:///{}'.format(self.eval_args['db'])
        for table in self.build_table_names():
            yield luigi.contrib.sqla.SQLAlchemyTarget(connection_string=conn_str,
                                                      target_table=table,
                                                      update_id='_'.join([self.task_id, table]))

    def requires(self):
        return self.clone(AlignTranscripts)

    def run(self):
        logger.info('Beginning Evaluate Transcript toil pipeline for {}.'.format(self.genome))
        toil_work_dir = os.path.join(self.work_dir, 'toil', 'transcript_evaluation', self.genome)
        toil_options = self.prepare_toil_options(toil_work_dir)
        results = classify(self.eval_args, toil_options)
        # results should be a dictionary of {table: dataframe}
        logger.info('Evaluate Transcript toil pipeline for {} completed'.format(self.genome))
        self.write_to_sql(results)


def parse_args():
    """
    If we are running from a local scheduler, we need to parse the arguments.
    :return: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--hal', required=True)
    parser.add_argument('--ref-genome', required=True)
    parser.add_argument('--annotation', required=True)
    parser.add_argument('--out-dir', default='./cat_output')
    parser.add_argument('--work-dir', default=tempfile.gettempdir())
    parser.add_argument('--target-genomes', nargs='+', default=None)
    # augustus TM(R) options
    parser.add_argument('--augustus', action='store_true')
    parser.add_argument('--augustus-species', default='human')
    parser.add_argument('--augustus-hints-db', default=None)
    parser.add_argument('--tm-cfg', default='augustus_cfgs/extrinsic.ETM1.cfg')
    parser.add_argument('--tmr-cfg', default='augustus_cfgs/extrinsic.ETM2.cfg')
    # augustus CGP options
    parser.add_argument('--augustus-cgp', action='store_true')
    parser.add_argument('--augustus-cgp-cfg', default=None)
    parser.add_argument('--augustus-cgp-param', default='augustus_cfgs/log_reg_parameters_default.cfg')
    parser.add_argument('--maf-chunksize', default=2500000, type=int)
    parser.add_argument('--maf-overlap', default=500000, type=int)
    parser.add_argument('--tm-to-hints-script', default='tools/transMap2hints.pl')
    # toil options
    parser.add_argument('--workDir', default=None)
    parser.add_argument('--batchSystem', default='singleMachine')
    parser.add_argument('--maxCores', default=16, type=int)
    parser.add_argument('--logLevel', default='WARNING')
    parser.add_argument('--cleanWorkDir', default='onSuccess')
    parser.add_argument('--parasolCommand', default=None)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    luigi.build([RunCat(**vars(args))], logging_conf_file='logging.cfg')
