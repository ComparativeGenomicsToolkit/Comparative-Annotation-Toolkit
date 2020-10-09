"""
Given a GFF3 file, print an uber-CDS BED file to stdout.
"""
import argparse
from tools.dataOps import flatten_list_of_lists
from tools.fileOps import TemporaryFilePath
from tools.procOps import run_proc
from tools.gff3 import convert_gff3_cmd
from tools.transcripts import get_gene_pred_dict, intervals_to_bed, group_transcripts_by_name2, Transcript
from tools.intervals import gap_merge_intervals


def parse_args():
    parser = argparse.ArgumentParser(usage="Given a GFF3 file, print an uber-CDS BED file to stdout.")
    parser.add_argument("gff", help="GFF file to parse")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    cmd = convert_gff3_cmd("/dev/null", args.gff)
    with TemporaryFilePath() as tmp:
        run_proc(cmd, stdout=tmp)
        tx_dict = get_gene_pred_dict(tmp)
    tx_dict = group_transcripts_by_name2(tx_dict.values())
    for gene, txs in tx_dict.items():
        cds_txs = [
            Transcript(tx.get_bed(new_start=tx.thick_start, new_stop=tx.thick_stop)) for tx in txs if tx.cds_size > 0
        ]
        if len(cds_txs) == 0:
            continue
        intervals = flatten_list_of_lists([x.exon_intervals for x in cds_txs])
        merged = gap_merge_intervals(intervals, 0)
        final_tx = intervals_to_bed(merged, name=gene)
        print("\t".join(final_tx.get_bed()))
