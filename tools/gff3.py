import collections
import pandas as pd
from . import fileOps
from . import transcripts
from . import misc


reserved_keys = ['gene_biotype',
                 'transcript_biotype',
                 'gene_type',
                 'transcript_type',
                 'gene_name',
                 'gene_id',
                 'transcript_id',
                 'transcript_name',
                 'ID',
                 'Name',
                 'Parent']


def parse_gff3(annotation_attrs, annotation_gp, is_external_reference=False):
    def parse_attrs(attrs):
        r = collections.defaultdict(dict)
        for tx_id, key, value in fileOps.iter_lines(attrs):
            r[tx_id][key] = value
        return r

    attrs_dict = parse_attrs(annotation_attrs)
    tx_dict = transcripts.get_gene_pred_dict(annotation_gp)
    tx_name_map = {x: y.name2 for x, y in tx_dict.items()}
    results = []
    for tx_id, gene_id in tx_name_map.items():
        d = attrs_dict[tx_id]
        gene_biotype = d.get('gene_biotype', d.get('gene_type', None))
        if gene_biotype is None:
            raise Exception("Did not find a gene biotype or gene type for {} (attrs={})".format(gene_id, d))
        tx_biotype = d.get('transcript_biotype', d.get('transcript_type', None))
        if tx_biotype is None:
            raise Exception("Did not find a transcript biotype or type for {} (attrs={})".format(tx_id, d))
        gene_name = d['gene_name']
        gene_id = d['gene_id']
        tx_id = d['transcript_id']
        tx_name = d['transcript_name']
        extra_tags = ';'.join(['{}={}'.format(x, y.replace(';', '%3B').replace('=', '%3D'))
                               for x, y in d.items() if x not in reserved_keys])
        if len(extra_tags) > 0:
            try:
                misc.parse_gff_attr_line(extra_tags)
            except:
                raise Exception(f'Error parsing extra tags in input GFF3 {extra_tags}')
        if is_external_reference is True:
            # hack to fix names
            gene_id = f'exRef-{gene_id}'
            tx_id = f'exRef-{tx_id}'
        results.append([gene_id, tx_id, tx_name, gene_name, gene_biotype, tx_biotype, extra_tags])
    df = pd.DataFrame(results, columns=['GeneId', 'TranscriptId', 'TranscriptName', 'GeneName',
                                        'GeneBiotype', 'TranscriptBiotype', 'ExtraTags'])
    df = df.set_index('TranscriptId')
    return df


def convert_gff3_cmd(annotation_attrs, annotation):
    cmd = ['gff3ToGenePred', '-rnaNameAttr=transcript_id', '-geneNameAttr=gene_id', '-honorStartStopCodons',
            '-refseqHacks',
           '-attrsOut={}'.format(annotation_attrs), annotation, '/dev/stdout']
    return cmd
