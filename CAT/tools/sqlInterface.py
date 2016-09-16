"""
Functions to interface with the sqlite databases produced by various steps of the annotation pipeline
"""
import collections
import sqlalchemy
import pandas as pd


###
# Useful functions for querying the databases
###


def read_attrs(db_path, table='annotation', index_col='TranscriptId'):
    """
    Read the attributes database file into a pandas DataFrame
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: pandas DataFrame
    """
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_path))
    return pd.read_sql_table(table, engine, index_col=index_col)


def get_transcript_gene_map(db_path, table='annotation', index_col='TranscriptId'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping transcript IDs to gene IDs.
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {tx_id: GeneId}
    """
    df = read_attrs(db_path, table, index_col)
    return dict(zip(df.index, df.GeneId))


def get_gene_transcript_map(db_path, table='annotation', index_col='TranscriptId'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping transcript IDs to gene IDs.
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {GeneId: [tx_id1, tx_id2, etc]}
    """
    df = read_attrs(db_path, table, index_col)
    r = collections.defaultdict(list)
    for tx_id, row in df.iterrows():
        r[row.GeneId].append(tx_id)
    return dict(r)  # don't return a defaultdict to prevent bugs


def get_transcript_biotype_map(db_path, table='annotation', index_col='TranscriptId'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping transcript IDs to their biotype
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {tx_id: tx_biotype}
    """
    df = read_attrs(db_path, table, index_col)
    return dict(zip(df.index, df.TranscriptBiotype))


def get_gene_biotype_map(db_path, table='annotation', index_col='TranscriptId'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping gene IDs to their biotype
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {tx_id: tx_biotype}
    """
    df = read_attrs(db_path, table, index_col)
    return dict(zip(df.GeneId, df.GeneBiotype))


def load_reference(ref_db_path):
    """
    Load the reference annotation table
    :param ref_db_path: path to reference genome database. Must have table 'annotation'
    :return: DataFrame
    """
    engine = sqlalchemy.create_engine('sqlite:///' + ref_db_path)
    df = pd.read_sql('annotation', engine, index_col=['GeneId', 'TranscriptId'])
    return df


def load_alignment_evaluation(db_path):
    """
    Loads the transMap alignment evaluation table
    :param db_path: path to genome database
    :return: DataFrame
    """
    engine = sqlalchemy.create_engine('sqlite:///' + db_path)
    df = pd.read_sql('alignment', engine)
    df['value'] = df['value'].apply(pd.to_numeric)
    return pd.pivot_table(df, index='AlignmentId', columns='classifier', values='value', fill_value=0)


def load_classifications(db_path, alignment_mode, transcript_modes):
    """
    Load all of the evaluation and metrics tables into one joined DataFrame
    :param db_path: path to genome database
    :param alignment_mode: one of ('CDS', 'mRNA')
    :param transcript_modes: List that contains modes ['transMap', 'augTM', 'augTMR', 'augCGP']
    :return: DataFrame
    """
    def load_evaluation(tx_type):
        """load evaluation table. Makes use of count() and group by to get the # of times the classifier failed"""
        table = '_'.join([alignment_mode, tx_type, 'Evaluation'])
        query = 'SELECT AlignmentId,TranscriptId,classifier,COUNT(*) FROM {} ' \
                'group by AlignmentId,TranscriptId,classifier'.format(table)
        df = pd.read_sql(query, engine)
        df.columns = ['AlignmentId', 'TranscriptId', 'classifier', 'value']
        return df

    def load_metrics(tx_type):
        """load metrics table"""
        table = '_'.join([alignment_mode, tx_type, 'Metrics'])
        query = 'SELECT AlignmentId,TranscriptId,classifier,value FROM {}'.format(table)
        df = pd.read_sql(query, engine)
        return df

    # we did not perform mRNA alignments on CGP
    if alignment_mode == 'mRNA':
        transcript_modes = list(set(transcript_modes) - {'augCGP'})
    engine = sqlalchemy.create_engine('sqlite:///' + db_path)
    dfs = [load_metrics(tx_mode) for tx_mode in transcript_modes]
    dfs.extend([load_evaluation(tx_mode) for tx_mode in transcript_modes])
    df = pd.concat(dfs, join='outer')
    df['value'] = df['value'].apply(pd.to_numeric)
    eval_df = pd.pivot_table(df, index=['AlignmentId', 'TranscriptId'], columns='classifier', values='value',
                             fill_value=0)
    return eval_df


def load_intron_vector(db_path, transcript_modes):
    """
    load the homGeneMapping intron vector, collapsing to a single value
    TODO: Make use of the deeper information present in the intron vector
    :param db_path: path to genome database
    :param transcript_modes: List that contains modes ['transMap', 'augTM', 'augTMR', 'augCGP']
    :return: DataFrame
    """
    def reduce_intron_vector(s):
        """intron vector is stored as a comma separated string. Apply this function to reduce to a integer count"""
        s.IntronVector = len([x for x in s.IntronVector.split(',') if x > 0])
        return s
    engine = sqlalchemy.create_engine('sqlite:///' + db_path)
    dfs = [pd.read_sql('_'.join([tx_mode, 'HgmIntronVector']), engine) for tx_mode in transcript_modes]
    dfs = [df.apply(reduce_intron_vector, axis=1) for df in dfs]
    df = pd.concat(dfs, join='outer')
    df.columns = ['AlignmentId', 'NumSupportedIntrons']
    return df
