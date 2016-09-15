"""
Functions to interface with the sqlite databases produced by various steps of the annotation pipeline
"""
import collections
import sqlalchemy
import pandas as pd


def read_attrs(db_path, table='annotation', index_col='tx_id'):
    """
    Read the attributes database file into a pandas DataFrame
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: pandas DataFrame
    """
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_path))
    return pd.read_sql_table(table, engine, index_col=index_col)


def get_transcript_gene_map(db_path, table='annotation', index_col='tx_id'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping transcript IDs to gene IDs.
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {tx_id: gene_id}
    """
    df = read_attrs(db_path, table, index_col)
    return dict(zip(df.index, df.gene_id))


def get_gene_transcript_map(db_path, table='annotation', index_col='tx_id'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping transcript IDs to gene IDs.
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {gene_id: [tx_id1, tx_id2, etc]}
    """
    df = read_attrs(db_path, table, index_col)
    r = collections.defaultdict(list)
    for tx_id, row in df.iterrows():
        r[row.gene_id].append(tx_id)
    return dict(r)  # don't return a defaultdict to prevent bugs


def get_transcript_biotype_map(db_path, table='annotation', index_col='tx_id'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping transcript IDs to their biotype
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {tx_id: tx_biotype}
    """
    df = read_attrs(db_path, table, index_col)
    return dict(zip(df.index, df.tx_biotype))


def get_gene_biotype_map(db_path, table='annotation', index_col='tx_id'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping gene IDs to their biotype
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {tx_id: tx_biotype}
    """
    df = read_attrs(db_path, table, index_col)
    return dict(zip(df.gene_id, df.gene_biotype))
