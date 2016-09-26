"""
Functions to interface with the sqlite databases produced by various steps of the annotation pipeline
"""
import collections

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
    GeneName = Column(Text)
    GeneBiotype = Column(Text)
    StartCodon = Column(Boolean)
    StopCodon = Column(Boolean)


class EvaluationColumns(object):
    """Mixin class for all TranscriptEvaluation module tables"""
    GeneId = Column(Text, primary_key=True)
    TranscriptId = Column(Text, primary_key=True)
    AlignmentId = Column(Text, primary_key=True)
    classifier = Column(Text, primary_key=True)
    chromosome = Column(Text)
    start = Column(Integer)
    stop = Column(Integer)
    strand = Column(Text)


class TmEval(EvaluationColumns, Base):
    """Table for evaluations from TransMapEvaluation module"""
    __tablename__ = 'TransMapEvaluation'


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


class CdsAugCgpEval(EvaluationColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusCGP"""
    __tablename__ = 'CDS_augCGP_Evaluation'


class MetricsColumns(object):
    """Mixin class for all TranscriptMetrics module tables"""
    GeneId = Column(Text, primary_key=True)
    TranscriptId = Column(Text, primary_key=True)
    AlignmentId = Column(Text, primary_key=True)
    classifier = Column(Text, primary_key=True)
    value = Column(Float)


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


class CdsAugCgpMetrics(MetricsColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusCGP"""
    __tablename__ = 'CDS_augCGP_Metrics'


class HgmColumns(object):
    """Mixin class for all homGeneMapping tables"""
    GeneId = Column(Text, primary_key=True)
    TranscriptId = Column(Text, primary_key=True)
    AlignmentId = Column(Text, primary_key=True)
    IntronVector = Column(Text)


class TmIntronSupport(HgmColumns, Base):
    """Table for intron support of transMap transcripts from homGeneMapping"""
    __tablename__ = 'transMap_HgmIntronVector'


class AugTmIntronSupport(HgmColumns, Base):
    """Table for intron support of AugustusTM transcripts from homGeneMapping"""
    __tablename__ = 'augTM_HgmIntronVector'


class AugTmrIntronSupport(HgmColumns, Base):
    """Table for intron support of AugustusTMR transcripts from homGeneMapping"""
    __tablename__ = 'augTMR_HgmIntronVector'


class AugCgpIntronSupport(HgmColumns, Base):
    """Table for intron support of AugustusCGP transcripts from homGeneMapping"""
    __tablename__ = 'augCGP_HgmIntronVector'


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
                  'augTMR': AugTmrIntronSupport, 'transMap': TmIntronSupport},
          'CDS': {'augCGP': {'metrics': CdsAugCgpMetrics, 'evaluation': CdsAugCgpEval},
                  'augTM': {'metrics': CdsAugTmMetrics, 'evaluation': CdsAugTmEval},
                  'augTMR': {'metrics': CdsAugTmrMetrics, 'evaluation': CdsAugTmrEval},
                  'transMap': {'metrics': CdsTmMetrics, 'evaluation': CdsTmEval}},
          'mRNA': {'augTM': {'metrics': MrnaAugTmMetrics, 'evaluation': MrnaAugTmEval},
                   'augTMR': {'metrics': MrnaAugTmrMetrics, 'evaluation': MrnaAugTmrEval},
                   'transMap': {'metrics': MrnaTmMetrics, 'evaluation': MrnaTmEval}}}


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
    df = read_attrs(db_path, table, index_col)
    r = collections.defaultdict(set)
    for tx_id, row in df.iterrows():
        r[row.GeneId].add(tx_id)
    return dict(r)  # don't return a defaultdict to prevent bugs


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
    df = pd.pivot_table(df, index=['TranscriptId', 'AlignmentId'], columns='classifier', values='value', fill_value=0)
    return df.reset_index()


###
# Load subsets of tables using the ORM
###


def _generic_gene_query(table, session, gene_id):
    """Generic function for selecting members of a table by gene_id"""
    query = session.query(table).filter(table.GeneId == gene_id)
    return pd.read_sql(query.statement, session.bind)


def _generic_transcript_query(table, session, tx_id):
    """Generic function for selecting members of a table by tx_id"""
    query = session.query(table).filter(table.TranscriptId == tx_id)
    return pd.read_sql(query.statement, session.bind)


def _generic_alignment_query(table, session, aln_id):
    """Generic function for selecting members of a table by aln_id"""
    query = session.query(table).filter(table.AlignmentId == aln_id)
    return pd.read_sql(query.statement, session.bind)


def load_evaluation(table, session, gene_id):
    """
    load evaluation entries for this gene. Makes use of count() and group by to get the # of times the classifier failed
    :param table: One of the evaluation tables
    :param session: Active sqlalchemy session.
    :param gene_id: Gene to query
    :return: DataFrame
    """
    assert any(table == cls for cls in (MrnaAugTmrEval, MrnaAugTmEval, MrnaTmEval,
                                        CdsAugCgpEval, CdsAugTmrEval, CdsAugTmEval, CdsTmEval))
    query = session.query(table.GeneId, table.TranscriptId, table.AlignmentId, table.classifier,
                          func.count(table.classifier).label('value')). \
        group_by(table.AlignmentId, table.TranscriptId, table.classifier). \
        filter(table.GeneId == gene_id)
    return pd.read_sql(query.statement, session.bind)


def load_metrics(table, session, gene_id):
    """
    load metrics entries for this gene. Wrapper for generic_gene_query.
    :param table: One of the metrics tables
    :param session: Active sqlalchemy session.
    :param gene_id: Gene to query
    :return: DataFrame
    """
    assert any(table == cls for cls in (MrnaAugTmrMetrics, MrnaAugTmMetrics, MrnaTmMetrics,
                                        CdsAugCgpMetrics, CdsAugTmrMetrics, CdsAugTmMetrics, CdsTmMetrics))
    return _generic_gene_query(table, session, gene_id)


def load_intron_vector(table, session, gene_id):
    """
    load intron vector entries for this gene. Wrapper for generic_gene_query.
    :param table: One of the intron vector tables
    :param session: Active sqlalchemy session.
    :param gene_id: Gene to query
    :return: DataFrame
    """
    assert any(table == cls for cls in (TmIntronSupport, AugCgpIntronSupport, AugTmIntronSupport,
                                        AugTmrIntronSupport))
    return _generic_gene_query(table, session, gene_id)
