"""
Ensembl GFF3 for non-mouse/human do not follow the format convention of having only gene/transcript/exon be the primary
feature type (column 2). This script fixes this and makes the GFF3 parseable by CAT.
"""
import argparse


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
            if 'ID=gene:' in l[-1]:
                l[2] = 'gene'
            elif 'ID=transcript:' in l[-1]:
                l[2] = 'transcript'
            args.output_gff3.write('\t'.join(l))


if __name__ == '__main__':
    main()
