"""
This module interfaces with the hints database produced for Augustus, providing a SQLAlchemy ORM access to it.
"""
import sqlalchemy
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import sessionmaker


def reflect_hints_db(db_path):
    """
    Reflect the database schema of the hints database, automapping the existing tables
    :param db_path: path to hints sqlite database
    :return: sqlalchemy.MetaData object, sqlalchemy.orm.Session object
    """
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_path))
    metadata = sqlalchemy.MetaData()
    metadata.reflect(bind=engine)
    Base = automap_base(metadata=metadata)
    Base.prepare()
    speciesnames = Base.classes.speciesnames
    seqnames = Base.classes.seqnames
    hints = Base.classes.hints
    featuretypes = Base.classes.featuretypes
    Session = sessionmaker(bind=engine)
    session = Session()
    return speciesnames, seqnames, hints, featuretypes, session


def get_rnaseq_hints(genome, chromosome, start, stop, speciesnames, seqnames, hints, featuretypes, session):
    """
    Extracts RNAseq hints from RNAseq hints database.
    :param genome: genome (table) to query
    :param chromosome: Chromosome to extract information from
    :param start: start position on chromosome
    :param stop: stop position in chromosome
    :param speciesnames: speciesnames Table from reflect_hints_db
    :param seqnames: seqnames Table from reflect_hints_db
    :param hints: hints Table from reflect_hints_db
    :param featuretypes: featuretypes Table from reflect_hints_db
    :param session: Session object from reflect_hints_db
    :return: GFF formatted string.
    """
    speciesid = session.query(speciesnames.speciesid).filter_by(speciesname=genome)
    seqnr = session.query(seqnames.seqnr).filter(
        sqlalchemy.and_(
            seqnames.speciesid.in_(speciesid),
            (seqnames.seqname == chromosome)))
    query = session.query(hints, featuretypes).filter(
            sqlalchemy.and_(
                hints.speciesid.in_(speciesid),
                hints.seqnr.in_(seqnr),
                hints.start >= start,
                hints.end <= stop,
                featuretypes.typeid == hints.type))
    hints = []
    for h, f in query:
        tags = 'pri=3;src={};mult={}'.format(h.esource, h.mult)
        l = [chromosome, h.source, f.typename, h.start + 1, h.end + 1, h.score, '.', '.', tags]
        hints.append('\t'.join(map(str, l)) + '\n')
    return '\n'.join(hints)
