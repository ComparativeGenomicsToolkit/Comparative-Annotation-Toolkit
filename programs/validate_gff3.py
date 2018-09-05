"""
This script allows you to test your GFF3 for compatibility with CAT before running the pipeline.

"""

import argparse
import tools.gff3
import tools.procOps
import tools.fileOps
import tools.transcripts
import collections

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('gff3', help='GFF3 to validate')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    with tools.fileOps.TemporaryFilePath() as attrs, tools.fileOps.TemporaryFilePath() as gp:
        cmd = tools.gff3.convert_gff3_cmd(attrs, args.gff3)
        tools.procOps.run_proc(cmd, stdout=gp)
        df = tools.gff3.parse_gff3(attrs, gp)
        tx_dict = tools.transcripts.get_gene_pred_dict(gp)
    assert len(tx_dict) == len(df)
    assert tx_dict.viewkeys() == set(df.index)
    genes = {x.name2 for x in tx_dict.itervalues()}
    assert genes == set(df.GeneId)
    print 'Found {} transcripts and {} genes'.format(len(tx_dict), len(genes))
    tmp = df.groupby('GeneId').first()
    print 'Found the following gene biotypes: {}'.format(collections.Counter(tmp.GeneBiotype))
    print 'Found the following transcript biotypes: {}'.format(collections.Counter(df.TranscriptBiotype))
    print 'Some example database rows:'
    print df.head(20)
    print 'Some example gene rows:'
    print tmp.head(10)