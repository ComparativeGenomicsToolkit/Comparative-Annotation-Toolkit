"""
Functions to interface with the sqlite databases produced by various steps of the annotation pipeline
"""
import transcripts

import pandas as pd
from sqlalchemy import Column, Integer, Text, Float, Boolean, func, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

###
# Data model
###


Base = declarative_base()


class Annotation(Base):
    """Table for the annotation table. Only exists in ref_genome"""
    __tablename__ = 'annotation'
    GeneId = Column(Text, primary_key=True)
    TranscriptId = Column(Text, primary_key=True)
    TranscriptName = Column(Text)
    GeneName = Column(Text)
    GeneBiotype = Column(Text)
    TranscriptBiotype = Column(Text)


class Bed12(object):
    """General table description for storing BED12 features"""
    chromosome = Column(Text)
    start = Column(Integer)
    stop = Column(Integer)
    name = Column(Text)
    score = Column(Integer)
    strand = Column(Text)
    thickStart = Column(Integer)
    thickStop = Column(Integer)
    rgb = Column(Text)
    blockCount = Column(Integer)
    blockSizes = Column(Text)
    blockStarts = Column(Text)


class EvaluationColumns(Bed12):
    """Mixin class for all TranscriptEvaluation module tables. Represents a bed12 with a leading ID column"""
    AlignmentId = Column(Text, primary_key=True)


class MrnaTmEval(EvaluationColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from transMap"""
    __tablename__ = 'mRNA_transMap_Evaluation'


class MrnaAugTmEval(EvaluationColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from AugustusTM"""
    __tablename__ = 'mRNA_augTM_Evaluation'


class MrnaAugTmrEval(EvaluationColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from AugustusTMR"""
    __tablename__ = 'mRNA_augTMR_Evaluation'


class CdsTmEval(EvaluationColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from transMap"""
    __tablename__ = 'CDS_transMap_Evaluation'


class CdsAugTmEval(EvaluationColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusTM"""
    __tablename__ = 'CDS_augTM_Evaluation'


class CdsAugTmrEval(EvaluationColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusTMR"""
    __tablename__ = 'CDS_augTMR_Evaluation'


class MetricsColumns(object):
    """Mixin class for all TranscriptMetrics module tables"""
    AlignmentId = Column(Text, primary_key=True)
    classifier = Column(Text)
    value = Column(Float)


class TmEval(MetricsColumns, Base):
    """Table for evaluations from TransMapEvaluation module"""
    __tablename__ = 'TransMapEvaluation'
    TranscriptId = Column(Text, primary_key=True)
    GeneId = Column(Text, primary_key=True)


class TmFilterEval(MetricsColumns, Base):
    """Table for evaluations from FilterTransMap module. This table is stored in a stacked format for simplicity."""
    __tablename__ = 'TransMapFilterEvaluation'
    GeneId = Column(Text, primary_key=True)
    TranscriptId = Column(Text, primary_key=True)
    AlignmentId = Column(Text, primary_key=True)
    GeneAlternateContigs = Column(Text)
    GeneAlternateLoci = Column(Text)
    CollapsedGeneNames = Column(Text)
    CollapsedGeneIds = Column(Text)
    Paralogy = Column(Text)
    UnfilteredParalogy = Column(Text)


class TmMetrics(MetricsColumns, Base):
    """Table for evaluations from TransMapMetrics module"""
    __tablename__ = 'TransMapMetrics'


class MrnaTmMetrics(MetricsColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from transMap"""
    __tablename__ = 'mRNA_transMap_Metrics'


class MrnaAugTmMetrics(MetricsColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from AugustusTM"""
    __tablename__ = 'mRNA_augTM_Metrics'


class MrnaAugTmrMetrics(MetricsColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from AugustusTMR"""
    __tablename__ = 'mRNA_augTMR_Metrics'


class CdsTmMetrics(MetricsColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from transMap"""
    __tablename__ = 'CDS_transMap_Metrics'


class CdsAugTmMetrics(MetricsColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusTM"""
    __tablename__ = 'CDS_augTM_Metrics'


class CdsAugTmrMetrics(MetricsColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusTMR"""
    __tablename__ = 'CDS_augTMR_Metrics'


class HgmColumns(object):
    """Mixin class for all homGeneMapping tables"""
    GeneId = Column(Text, primary_key=True)
    TranscriptId = Column(Text, primary_key=True)
    AlignmentId = Column(Text, primary_key=True)
    AllSpeciesIntronRnaSupport = Column(Text)
    AllSpeciesExonRnaSupport = Column(Text)
    IntronRnaSupport = Column(Text)
    ExonRnaSupport = Column(Text)
    IntronAnnotSupport = Column(Text)
    CdsAnnotSupport = Column(Text)
    ExonAnnotSupport = Column(Text)


class TmIntronSupport(HgmColumns, Base):
    """Table for intron support of transMap transcripts from homGeneMapping"""
    __tablename__ = 'transMap_Hgm'


class AugTmIntronSupport(HgmColumns, Base):
    """Table for intron support of AugustusTM transcripts from homGeneMapping"""
    __tablename__ = 'augTM_Hgm'


class AugTmrIntronSupport(HgmColumns, Base):
    """Table for intron support of AugustusTMR transcripts from homGeneMapping"""
    __tablename__ = 'augTMR_Hgm'


class AugCgpIntronSupport(HgmColumns, Base):
    """Table for intron support of AugustusCGP transcripts from homGeneMapping"""
    __tablename__ = 'augCGP_Hgm'


class AugPbIntronSupport(HgmColumns, Base):
    """Table for intron support of AugustusPB transcripts from homGeneMapping"""
    __tablename__ = 'augPB_Hgm'


class ExRefIntronSupport(HgmColumns, Base):
    """Table for intron support of External reference transcripts from homGeneMapping"""
    __tablename__ = 'ExRef_Hgm'


class AlternativeGeneIdColumns(object):
    """mixin class for AlternativeGenes"""
    TranscriptId = Column(Text, primary_key=True)
    AssignedGeneId = Column(Text)
    AlternativeGeneIds = Column(Text)
    ResolutionMethod = Column(Text)


class AugCgpAlternativeGenes(AlternativeGeneIdColumns, Base):
    """Table for recording a list of alternative parental genes for CGP"""
    __tablename__ = 'augCGP_AlternativeGenes'


class AugPbAlternativeGenes(AlternativeGeneIdColumns, Base):
    """Table for recording a list of alternative parental genes for IsoSeq"""
    __tablename__ = 'augPB_AlternativeGenes'


class ExRefAlternativeGenes(AlternativeGeneIdColumns, Base):
    """Table for recording a list of alternative parental genes for external references"""
    __tablename__ = 'ExRef_AlternativeGenes'


class IsoSeqExonStructures(Bed12, Base):
    """Table for recording all distinct exon structures present in a IsoSeq hints file"""
    __tablename__ = 'IsoSeqExonStructures'
    index = Column(Integer, primary_key=True)


###
# Wrapper functions for setting up sessions
###


def start_session(db_path):
    """basic script for starting a session"""
    engine = create_engine('sqlite:///' + db_path)
    Session = sessionmaker(bind=engine)
    return Session()


###
# Dictionary mapping tables to their respective transcript/alignment modes
###


tables = {'hgm': {'augCGP': AugCgpIntronSupport, 'augTM': AugTmIntronSupport,
                  'augTMR': AugTmrIntronSupport, 'transMap': TmIntronSupport,
                  'augPB': AugPbIntronSupport, 'exRef': ExRefIntronSupport},
          'CDS': {'augTM': {'metrics': CdsAugTmMetrics, 'evaluation': CdsAugTmEval},
                  'augTMR': {'metrics': CdsAugTmrMetrics, 'evaluation': CdsAugTmrEval},
                  'transMap': {'metrics': CdsTmMetrics, 'evaluation': CdsTmEval}},
          'mRNA': {'augTM': {'metrics': MrnaAugTmMetrics, 'evaluation': MrnaAugTmEval},
                   'augTMR': {'metrics': MrnaAugTmrMetrics, 'evaluation': MrnaAugTmrEval},
                   'transMap': {'metrics': MrnaTmMetrics, 'evaluation': MrnaTmEval}},
          'alt_names': {'exRef': ExRefAlternativeGenes,
                        'augPB': AugPbAlternativeGenes,
                        'augCGP': AugCgpAlternativeGenes}}


###
# Attributes functions -- read data from the annotation table
###


def read_attrs(db_path, table=Annotation.__tablename__, index_col='TranscriptId'):
    """
    Read the attributes database file into a pandas DataFrame
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: pandas DataFrame
    """
    engine = create_engine('sqlite:///{}'.format(db_path))
    return pd.read_sql_table(table, engine, index_col=index_col)


def get_transcript_gene_map(db_path, table=Annotation.__tablename__, index_col='TranscriptId'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping transcript IDs to gene IDs.
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {tx_id: GeneId}
    """
    df = read_attrs(db_path, table, index_col)
    return dict(zip(df.index, df.GeneId))


def get_gene_transcript_map(db_path, table=Annotation.__tablename__, index_col='TranscriptId'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping transcript IDs to gene IDs.
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {GeneId: [tx_id1, tx_id2, etc]}
    """
    df = read_attrs(db_path, table, index_col).reset_index()
    r = {}
    for gene_id, s in df.groupby('GeneId'):
        r[gene_id] = s.TranscriptId.tolist()
    return r


def get_transcript_biotype_map(db_path, table=Annotation.__tablename__, index_col='TranscriptId'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping transcript IDs to their biotype
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {tx_id: tx_biotype}
    """
    df = read_attrs(db_path, table, index_col)
    return dict(zip(df.index, df.TranscriptBiotype))


def get_gene_biotype_map(db_path, table=Annotation.__tablename__, index_col='TranscriptId'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping gene IDs to their biotype
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {tx_id: tx_biotype}
    """
    df = read_attrs(db_path, table, index_col)
    return dict(zip(df.GeneId, df.GeneBiotype))


def get_transcript_biotypes(db_path, table=Annotation):
    """
    Returns a set of transcript biotypes seen in this annotation set
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :return: dictionary {tx_id: tx_biotype}
    """
    session = start_session(db_path)
    query = session.query(table.TranscriptBiotype).distinct()
    return {x[0] for x in query.all()}


def get_gene_biotypes(db_path, table=Annotation):
    """
    Returns a set of transcript biotypes seen in this annotation set
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :return: dictionary {tx_id: tx_biotype}
    """
    session = start_session(db_path)
    query = session.query(table.GeneBiotype).distinct()
    return {x[0] for x in query.all()}


###
# Loading entire tables
###


def load_annotation(ref_db_path):
    """
    Load the reference annotation table
    :param ref_db_path: path to reference genome database. Must have table Annotation.__tablename__
    :return: DataFrame
    """
    engine = create_engine('sqlite:///' + ref_db_path)
    df = pd.read_sql_table(Annotation.__tablename__, engine)
    return df


def load_alignment_evaluation(db_path):
    """
    Loads the transMap alignment evaluation table
    :param db_path: path to genome database
    :return: DataFrame
    """
    engine = create_engine('sqlite:///' + db_path)
    df = pd.read_sql_table(TmEval.__tablename__, engine)
    df = pd.pivot_table(df, index=['TranscriptId', 'AlignmentId'], columns='classifier', values='value')
    return df.reset_index()


def load_filter_evaluation(db_path):
    """
    Loads the transMap alignment filtering evaluation table
    :param db_path: path to genome database
    :return: DataFrame
    """
    engine = create_engine('sqlite:///' + db_path)
    return pd.read_sql_table(TmFilterEval.__tablename__, engine)


def load_isoseq_txs(db_path):
    """
    Loads the table IsoSeqExonStructures, constructing actual ChromosomeInterval objects.
    :param db_path: path to genome db
    :return: list of Transcript objects
    """
    engine = create_engine('sqlite:///' + db_path)
    df = pd.read_sql_table(IsoSeqExonStructures.__tablename__, engine, index_col='index')
    txs = [transcripts.Transcript(list(s)) for _, s in df.iterrows()]
    return txs


def load_evaluation(table, session):
    """
    load evaluation entries for this gene. Makes use of count() and group by to get the # of times the classifier failed
    :param table: One of the evaluation tables
    :param session: Active sqlalchemy session.
    :return: DataFrame
    """
    assert any(table == cls for cls in (MrnaAugTmrEval, MrnaAugTmEval, MrnaTmEval,
                                        CdsAugTmrEval, CdsAugTmEval, CdsTmEval))
    query = session.query(table.AlignmentId, table.name, func.count(table.name).label('value')). \
        group_by(table.AlignmentId, table.name)
    return pd.read_sql(query.statement, session.bind)


def load_metrics(table, session):
    """
    load metrics entries for this gene. Wrapper for generic_gene_query.
    :param table: One of the metrics tables
    :param session: Active sqlalchemy session.
    :return: DataFrame
    """
    assert any(table == cls for cls in (MrnaAugTmrMetrics, MrnaAugTmMetrics, MrnaTmMetrics,
                                        CdsAugTmrMetrics, CdsAugTmMetrics, CdsTmMetrics))
    query = session.query(table)
    return pd.read_sql(query.statement, session.bind)


def load_intron_vector(table, session):
    """
    load intron vector entries for this gene. Wrapper for generic_gene_query.
    :param table: One of the intron vector tables
    :param session: Active sqlalchemy session.
    :return: DataFrame
    """
    assert any(table == cls for cls in (TmIntronSupport, AugCgpIntronSupport, AugTmIntronSupport, AugPbIntronSupport,
                                        AugTmrIntronSupport))
    query = session.query(table)
    return pd.read_sql(query.statement, session.bind)


def load_alternatives(table, session):
    """
    load AugustusCGP/PB parental assignment + alternative parents
    :param table: Either AugCgpAlternativeGenes or AugPbAlternativeGenes
    :param session: Active sqlalchemy session.
    :return: DataFrame
    """
    assert table == AugCgpAlternativeGenes or table == AugPbAlternativeGenes
    query = session.query(table)
    return pd.read_sql(query.statement, session.bind)


###
# Stats functions
###

def load_luigi_stats(db_path, table):
    """
    Loads the luigi stats from the stats db
    :param db_path: path to database
    :return: DataFrame
    """
    engine = create_engine('sqlite:///' + db_path)
    return pd.read_sql_table(table, engine)
