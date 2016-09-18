"""
Functions to interface with the sqlite databases produced by various steps of the annotation pipeline
"""
import collections
import sqlalchemy
import pandas as pd


###
# Attributes functions
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


###
# Loading classification tables
###


def load_reference(ref_db_path):
    """
    Load the reference annotation table
    :param ref_db_path: path to reference genome database. Must have table 'annotation'
    :return: DataFrame
    """
    engine = sqlalchemy.create_engine('sqlite:///' + ref_db_path)
    df = pd.read_sql('annotation', engine)
    return df


def load_alignment_evaluation(db_path):
    """
    Loads the transMap alignment evaluation table
    :param db_path: path to genome database
    :return: DataFrame
    """
    engine = sqlalchemy.create_engine('sqlite:///' + db_path)
    df = pd.read_sql('TransMapEvaluation', engine)
    df['value'] = df['value'].apply(pd.to_numeric)
    df = pd.pivot_table(df, index=['TranscriptId', 'AlignmentId'], columns='classifier', values='value', fill_value=0)
    return df.reset_index()


###
# These functions require external information to create their dataframes
###


def load_classifications(db_path, alignment_mode, transcript_modes, ref_tx_dict):
    """
    Load all of the evaluation and metrics tables into one joined DataFrame
    :param db_path: path to genome database
    :param alignment_mode: one of ('CDS', 'mRNA')
    :param transcript_modes: List that contains modes ['transMap', 'augTM', 'augTMR', 'augCGP']
    :param ref_tx_dict: dictionary of GenePredTranscript objects representing the reference transcripts
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

    def add_intron_exon_counts():
        """based on alignment mode, produce a DataFrame of the number of reference introns/exons"""
        r = []
        for ref_tx_id, ref_tx in ref_tx_dict.iteritems():
            if alignment_mode == 'mRNA':
                r.append([ref_tx_id, len(ref_tx.exon_intervals), len(ref_tx.intron_intervals)])
            else:
                r.append([ref_tx_id, ref_tx.num_coding_exons, ref_tx.num_coding_introns])
        df = pd.DataFrame(r)
        df.columns = ['TranscriptId', 'NumReferenceExons', 'NumReferenceIntrons']
        return df

    # we did not perform mRNA alignments on CGP
    if alignment_mode == 'mRNA':
        transcript_modes = list(set(transcript_modes) - {'augCGP'})
    engine = sqlalchemy.create_engine('sqlite:///' + db_path)
    dfs = [load_metrics(tx_mode) for tx_mode in transcript_modes]
    dfs.extend([load_evaluation(tx_mode) for tx_mode in transcript_modes])
    df = pd.concat(dfs, join='outer')
    # we have to convert the value column to numeric for pivot to work
    df['value'] = df['value'].apply(pd.to_numeric)
    eval_df = pd.pivot_table(df, index=['AlignmentId', 'TranscriptId'], columns='classifier', values='value',
                             fill_value=0)
    # bring in the intron counts
    intron_df = add_intron_exon_counts()
    return pd.merge(eval_df.reset_index(), intron_df, on='TranscriptId')


def load_intron_vector(db_path, aln_mode, transcript_modes, tx_dict):
    """
    load the homGeneMapping intron vector, collapsing to a single value.
    We pass aln_mode/tx to only count CDS introns if we are in CDS mode
    TODO: Make use of the deeper information present in the intron vector
    :param db_path: path to genome database
    :param aln_mode: One of ('CDS', 'mRNA')
    :param transcript_modes: List that contains modes ['transMap', 'augTM', 'augTMR', 'augCGP']
    :param tx_dict: dictionary of GenePredTranscript objects representing the target transcripts
    :return: DataFrame
    """
    def reduce_intron_vector(s, aln_mode, tx, aln_id):
        """intron vector is stored as a comma separated string. Reduce this, taking aln_mode into account"""
        num_supported = 0
        scores = map(int, list(s)[0].split(','))
        for intron, score in zip(*[tx.intron_intervals, scores]):
            if aln_mode == 'CDS' and not intron.subset(tx.coding_interval):  # don't look at this intron
                continue
            if score == 0:  # this intron is not supported
                continue
            num_supported += 1
        return aln_id, num_supported

    engine = sqlalchemy.create_engine('sqlite:///' + db_path)
    dfs = [pd.read_sql('_'.join([tx_mode, 'HgmIntronVector']), engine) for tx_mode in transcript_modes]
    df = pd.concat(dfs, join='outer')
    df = df.set_index('AlignmentId')
    r = [reduce_intron_vector(s, aln_mode, tx_dict[aln_id], aln_id) for aln_id, s in df.iterrows()]
    df = pd.DataFrame(r)
    df.columns = ['AlignmentId', 'NumSupportedIntrons']
    return df
