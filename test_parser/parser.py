# test for gff parser
from pandas.util.testing import assert_frame_equal
import pandas as pd
import sys

def run(annotation_gp,annotation_attrs):
    df = pd.read_csv(annotation_attrs, sep='\t', names=['transcript_id', 'key', 'value'], header=None)
    results = []
    if 'gbkey' in set(df.key):
        for tx_id, d in df.groupby('transcript_id'):
            d = dict(zip(d.key, d.value))
            if d['gbkey'] == 'Gene':
	        continue
            elif d['gbkey'] in ['tRNA', 'rRNA']:
                gene_biotype = tx_biotype = d['gbkey']
            elif d['gbkey'] == 'mRNA':
                # hacky check because of lack of biotype features on transcript-level feature
                if 'pseudo' in d and d['pseudo'] == 'true':
                    gene_biotype = tx_biotype = 'pseudogene'
                else:
                    gene_biotype = tx_biotype = 'protein_coding'
            else:
                gene_biotype = tx_biotype = d['gbkey']
            gene_name = d.get('gene', 'Name')
            gene_id = d.get('Parent', tx_id)
            tx_name = d.get('product', tx_id)
            results.append([gene_id, tx_id, tx_name, gene_name, gene_biotype, tx_biotype])
    else:  # this is not a NCBI GFF3
        for tx_id, d in df.groupby('transcript_id'):
            d = dict(zip(d.key, d.value))
            if 'transcript_id' not in d:
                continue
            else:
                if 'biotype' in d:  # possibly Ensembl
                    gene_biotype = tx_biotype = d['biotype']
                elif 'gene_type' in d:  # probably Gencode
                    gene_biotype = d['gene_type']
                    tx_biotype = d['transcript_type']
                else:
                    raise InvalidInputException('Could not parse biotype for {}. Values: {}'.format(tx_id, d))
                # Ensembl formats their GFF3 with the format ID=transcript:XXX, while Gencode doesn't have the
                # extraneous transcript: portion.
                # Gencode also includes the gene name on the transcript level, so it is carried over.
                # Ensembl does not do this, but we can infer this via the regular schema Name-Version
                # However, Ensembl also does not always include a Name tag, so we have to account for this as well
                if 'transcript' in d['ID']:  # probably Ensembl
                    gene_id = d['Parent'].replace('gene:', '')
                    if 'Name' in d:
                        gene_name = d['Name'].split('-')[0]
                        tx_name = d['Name']
                    else:  # no names here, just use IDs
                        gene_name = gene_id
                        tx_name = tx_id
                elif 'gene_name' in d and 'gene_id' in d and 'transcript_name' in d:  # Gencode
                    gene_name = d['gene_name']
                    gene_id = d['gene_id']
                    tx_name = d['transcript_name']
                else:  # ambiguous type, hope for the best here
                    if 'gene' in d:
                        gene_name = d['gene']
                    elif 'Name' in d:
                        gene_name = d['Name']
                    else:
                        gene_name = d['Parent']
                    gene_id = d['Parent']
                    TX_name = d.get('product', tx_id)
            results.append([gene_id, tx_id, tx_name, gene_name, gene_biotype, tx_biotype])
    df = pd.DataFrame(results, columns=['GeneId', 'TranscriptId', 'TranscriptName', 'GeneName',
                                            'GeneBiotype', 'TranscriptBiotype'])
    df = df.set_index('TranscriptId')
    tot_genes = len(open(annotation_gp).readlines())
    if tot_genes != len(df):
        raise InvalidInputException('The number of genes parsed from the attrs file is not the same number as '
                                    'in the genePred. This is a parser failure. Contact Ian and make him fix it.')
    return df

def test_ncbi():
   test = run("ncbiGp","ncbiAttrs")
   ncbi = pd.read_csv("ncbiTest",sep='\t',index_col=0)
   assert_frame_equal(test,ncbi,check_dtype=False)

def test_gencode():
   test = run("gencodeGp","gencodeAttrs")
   gencode = pd.read_csv("gencodeTest",sep='\t',index_col=0)
   assert_frame_equal(test,gencode,check_dtype=False)

def test_ensembl():
   test = run("ensemblGp","ensemblAttrs")
   ensembl = pd.read_csv("ensemblTest",sep='\t',index_col=0)
   assert_frame_equal(test,ensembl,check_dtype=False)

