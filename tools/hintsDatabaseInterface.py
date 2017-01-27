"""
This module interfaces with the hints database produced for Augustus, providing a SQLAlchemy ORM access to it.
"""
import sqlalchemy
from sqlalchemy.pool import NullPool
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import sessionmaker


def reflect_hints_db(db_path):
    """
    Reflect the database schema of the hints database, automapping the existing tables

    The NullPool is used to avoid concurrency issues with luigi. Using this activates pooling, but since sqlite doesn't
    really support pooling, what effectively happens is just that it locks the database and the other connections wait.

    :param db_path: path to hints sqlite database
    :return: sqlalchemy.MetaData object, sqlalchemy.orm.Session object
    """
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_path), poolclass=NullPool)
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
        # add 1 to both start and end to shift to 1-based
        l = [chromosome, h.source, f.typename, h.start + 1, h.end + 1, h.score, '.', '.', tags]
        hints.append('\t'.join(map(str, l)) + '\n')
    return ''.join(hints)


def get_wiggle_hints(genome, speciesnames, seqnames, hints, session):
    """
    Extracts all wiggle hints for a genome to a BED format.
    :param genome: genome (table) to query
    :param speciesnames: speciesnames Table from reflect_hints_db
    :param seqnames: seqnames Table from reflect_hints_db
    :param hints: hints Table from reflect_hints_db
    :param session: Session object from reflect_hints_db
    :return: iterator of BED format lists
    """
    speciesid = session.query(speciesnames.speciesid).filter_by(speciesname=genome)
    seqs = {x.seqnr: x.seqname for x in session.query(seqnames).filter_by(speciesid=speciesid)}
    # chunk up the genome to reduce memory usage
    for seqnr, seqname in seqs.iteritems():
        query = session.query(hints.start, hints.end, hints.score).filter(
                sqlalchemy.and_(hints.speciesid.in_(speciesid), hints.source == 'w2h', hints.seqnr == seqnr))
        for start, end, score in query:
            # add 1 to end to convert to half-open interval
            yield seqname, start, end + 1, score


def hints_db_has_rnaseq(db_path, genome=None):
    """
    Determines if the hints DB has RNAseq. Is done by querying for one b2h or w2h in hints
    :param db_path: path to database
    :param genome: set this to query a specific genome instead of the database in general
    :return: boolean
    """
    speciesnames, seqnames, hints, featuretypes, session = reflect_hints_db(db_path)
    query = session.query(hints).filter(sqlalchemy.or_(hints.source == 'w2h', hints.source == 'b2h'))
    if genome is not None:
        speciesid = session.query(speciesnames.speciesid).filter_by(speciesname=genome)
        query = query.filter(hints.speciesid == speciesid)
    r = query.first() is not None
    session.close()
    return r


def genome_has_no_wiggle_hints(db_path, genome):
    """
    Determines if the hints db for a specific genome has wiggle hints
    :param db_path: path to database
    :param genome: genome in question
    :return: boolean
    """
    speciesnames, seqnames, hints, featuretypes, session = reflect_hints_db(db_path)
    query = session.query(hints).filter(hints.source == 'w2h')
    speciesid = session.query(speciesnames.speciesid).filter_by(speciesname=genome)
    query = query.filter(hints.speciesid == speciesid)
    r = query.first() is None
    session.close()
    return r


def hints_db_has_annotation(db_path, genome=None):
    """
    Determines if the hints DB has annotation. Is done by querying for a2h in hints
    :param db_path: path to database
    :param genome: set this to query a specific genome instead of the database in general
    :return: boolean
    """
    speciesnames, seqnames, hints, featuretypes, session = reflect_hints_db(db_path)
    query = session.query(hints).filter(hints.source == 'a2h')
    if genome is not None:
        speciesid = session.query(speciesnames.speciesid).filter_by(speciesname=genome)
        query = query.filter(hints.speciesid == speciesid)
    r = query.first() is not None
    session.close()
    return r
