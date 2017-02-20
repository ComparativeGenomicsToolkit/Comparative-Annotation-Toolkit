"""
NCBI GFF3 files have a few problems.

1) They use specific transcript level features such as mRNA.
2) They do not have the gene_id field.
3) They use gene_biotype instead of biotype.
4) They have no transcript biotypes.

"""
import argparse
from tools.misc import parse_gff_attr_line


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('gff3', type=argparse.FileType('r'))
    parser.add_argument('output_gff3', type=argparse.FileType('w'))
    return parser.parse_args()


def main():
    args = parse_args()
    for l in args.gff3:
        if l.startswith('#'):
            args.output_gff3.write(l)
        else:
            l = l.split('\t')
            if 'gbkey=Gene' in l[-1]:
                # fix gene level records
                attrs = parse_gff_attr_line(l[-1])
                # fix biotype
                attrs['biotype'] = attrs['gene_biotype']
                del attrs['gene_biotype']
                # fix gene id
                dbxref = attrs['Dbxref'].split(',')
                dbxref = dict(x.split(':') for x in dbxref)
                attrs['gene_id'] = dbxref['GeneID']
                del attrs['Dbxref']
                l[-1] = ';'.join(['='.join(i) for i in attrs.iteritems()]) + '\n'
                l[2] = 'gene'
            elif 'gbkey=mRNA' in l[-1] and l[2] not in ['CDS', 'exon', 'UTR']:
                l[2] = 'transcript'
            args.output_gff3.write('\t'.join(l))


if __name__ == '__main__':
    main()
