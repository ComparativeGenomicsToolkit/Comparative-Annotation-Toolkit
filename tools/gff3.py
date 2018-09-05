import collections
import pandas as pd
import fileOps
import transcripts

def parse_gff3(annotation_attrs, annotation_gp):
    def parse_attrs(attrs):
        r = collections.defaultdict(dict)
        for tx_id, key, value in fileOps.iter_lines(attrs):
            r[tx_id][key] = value
        return r

    attrs_dict = parse_attrs(annotation_attrs)
    tx_dict = transcripts.get_gene_pred_dict(annotation_gp)
    tx_name_map = {x: y.name2 for x, y in tx_dict.iteritems()}
    results = []
    for tx_id, gene_id in tx_name_map.iteritems():
        d = attrs_dict[tx_id]
        if 'gbkey' in d:  # NCBI
            if d['gbkey'] == 'mRNA':
                # hacky check because of lack of biotype features on transcript-level features
                if 'pseudo' in d and d['pseudo'] == 'true':
                    gene_biotype = tx_biotype = 'pseudogene'
                else:
                    gene_biotype = tx_biotype = 'protein_coding'
            elif d['gbkey'] == 'CDS':  # this is a transcript missing a transcript-level feature
                gene_biotype = tx_biotype = 'protein_coding'
            else:
                gene_biotype = tx_biotype = d['gbkey']
            if 'gene' in d:
                gene_name = d['gene']
            elif 'Name' in d:
                gene_name = d['Name']
            else:
                gene_name = d.get('Parent', 'ID')
            tx_name = d.get('product', tx_id)
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
                tx_name = d['transcript_name']
            else:  # ambiguous type, hope for the best here
                if 'gene' in d:
                    gene_name = d['gene']
                elif 'Name' in d:
                    gene_name = d['Name']
                else:
                    gene_name = d['Parent']
                tx_name = d.get('product', tx_id)
        results.append([gene_id, tx_id, tx_name, gene_name, gene_biotype, tx_biotype])
    df = pd.DataFrame(results, columns=['GeneId', 'TranscriptId', 'TranscriptName', 'GeneName',
                                        'GeneBiotype', 'TranscriptBiotype'])
    df = df.set_index('TranscriptId')
    return df


def convert_gff3_cmd(annotation_attrs, annotation):
    cmd = ['gff3ToGenePred', '-rnaNameAttr=transcript_id', '-geneNameAttr=gene_id', '-honorStartStopCodons',
           '-attrsOut={}'.format(annotation_attrs), annotation, '/dev/stdout']
    return cmd