#!/usr/bin/env python
"""Developed against Rnor6 RefSeq"""
import argparse
from BCBio import GFF
from collections import Counter, OrderedDict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_gff3', help='RefSeq GFF3 file')
    parser.add_argument('output_gff3', help='Output GFF3')
    return parser.parse_args()


def add_tags_to_feature(feature):
    if feature.type not in ['ncRNA', 'tRNA', 'rRNA', 'pseudogene', 'lnc_RNA',
                   'rRNA', 'mRNA', 'C_gene_segment',
                   'J_gene_segment', 'V_gene_segment',
                   'primary_transcript', 'miRNA', 'transcript', 'CDS', 'exon', 'gene']:
        return
    ids = dict([x.split(':') for x in feature.qualifiers.get('dbxref', feature.qualifiers['Dbxref'])])
    if feature.type == 'gene':
        feature.qualifiers['gene_id'] = [ids['GeneID']]
        return
    feature.qualifiers['gene_id'] = [ids['GeneID']]
    feature.qualifiers['gene_name'] = feature.qualifiers.get('gene', feature.qualifiers['gene_id'])
    tx_id = feature.qualifiers.get('transcript_id', feature.qualifiers['ID'])[0]
    if '-' in tx_id:
        tx_id = tx_id.split('-', 1)[1]
    feature.qualifiers['transcript_id'] = [tx_id]
    feature.qualifiers['transcript_name'] = feature.qualifiers.get('Name', feature.qualifiers['transcript_id'])
    if feature.type in ['mRNA', 'CDS']:
        gene_biotype = transcript_biotype = ['protein_coding']
    elif feature.type == 'primary_transcript':
        gene_biotype = ['miRNA']
        transcript_biotype = feature.type
    else:
        gene_biotype = transcript_biotype = [feature.type]
    feature.qualifiers['gene_biotype'] = gene_biotype
    feature.qualifiers['transcript_biotype'] = transcript_biotype


def construct_new_qualifiers(feature):
    new_qualifiers = OrderedDict()
    for key, val in feature.qualifiers.items():
        # no upper case keys unless it is ID or Parent or Name
        if key not in ['ID', 'Parent', 'Name']:
            key = key.lower()
        # collapse to a single item
        # replace all semicolons
        if len(val) > 1:
            val = [' '.join([x.replace(';', '%3B').replace('=', '%3D') for x in val])]
        new_qualifiers[key] = val
    # clean up and make parseable
    for key, val in new_qualifiers.items():
        if sum(len(x) for x in val) == 0:
            new_qualifiers[key] = 'True'
    return new_qualifiers


def feature_traversal(feature):
    yield feature
    for sub_feature in feature.sub_features:
        yield from feature_traversal(sub_feature)


if __name__ == '__main__':
    args = parse_args()
    records = list(GFF.parse(args.input_gff3))
    for seqrecord in records:
        for parent_feature in seqrecord.features:
            for feature in feature_traversal(parent_feature):
                try:
                    add_tags_to_feature(feature)
                except KeyError:
                    assert False, feature.qualifiers
                new_qualifiers = construct_new_qualifiers(feature)
                feature.qualifiers = new_qualifiers

    with open(args.output_gff3, 'w') as fh:
        GFF.write(records, fh)
