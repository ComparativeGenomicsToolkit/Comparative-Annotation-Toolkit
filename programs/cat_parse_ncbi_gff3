#!/usr/bin/env python
"""
Convert a eukaryotic NCBI GFF3 to be CAT compatible
"""
import gffutils
import argparse
from copy import deepcopy


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_gff3', help='NCBI GFF3 file')
    parser.add_argument('output_gff3', help='Output GFF3')
    return parser.parse_args()


def update_gene_attributes(gene):
    """
    Updates attributes for genes. Assumes that the gene is tagged by /gene, and may or may not have a gene_id.
    """
    gene.attributes["gene_name"] = gene.attributes["gene"]
    del gene.attributes["gene"]

    if "gene_id" not in gene.attributes:
        gene.attributes["gene_id"] = gene.attributes["gene_name"]
    assert "gene_biotype" in gene.attributes


def update_transcript_attributes(gene, transcript, i):
    """
    Updates attributes for transcripts. Uses information from the gene.

    Assumes that transcript_biotype == gene_biotype.

    Ensures that the transcript level feature has the required keys:

    1. gene_id
    2. gene_name
    3. transcript_id
    4. transcript_name
    5. transcript_biotype
    6. gene_biotype

    Args:
        gene: Current gene.
        transcript: Current transcript.
        i: Number of this transcript, used to create unique identifiers if necessary.
    """
    if "transcript_id" not in transcript.attributes:
        transcript.attributes["transcript_id"] = [f"transcript-{transcript.id.replace('id-', '')}-{i}"]

    # propagate gene information to transcript information
    transcript.attributes["transcript_name"] = gene.attributes["gene_name"]
    transcript.attributes["gene_id"] = gene.attributes["gene_id"]
    transcript.attributes["gene_name"] = gene.attributes["gene_name"]

    transcript.attributes["gene_biotype"] = gene.attributes["gene_biotype"]
    transcript.attributes["transcript_biotype"] = gene.attributes["gene_biotype"]


def infer_transcript(gene, i):
    """
    Creates a transcript from a gene. Sets up relationships.

    Args:
        gene: A feature.
        i: Index of this transcript.

    Returns:
        A new copy of the Feature converted to be a transcript.
    """
    tx = deepcopy(gene)
    tx.featuretype = "transcript"
    tx.attributes["Parent"] = gene.attributes["ID"]
    new_id = [f"transcript-{tx.attributes['ID'][0].replace('id-', '').replace('gene-', '')}-{i}"]
    tx.attributes["ID"] = tx.attributes["Name"] = tx.attributes["transcript_id"] = tx.attributes["transcript_name"] = new_id
    tx.attributes["transcript_biotype"] = tx.attributes["gene_biotype"]

    return tx


if __name__ == "__main__":
    args = parse_args()
    db = gffutils.create_db(args.input_gff3, dbfn=":memory:", merge_strategy="create_unique")
    with open(args.output_gff3, "w") as fh:
        print("##gff-version 3", file=fh)

        for gene in db.features_of_type(["gene", "pseudogene"]):

            update_gene_attributes(gene)
            print(gene, file=fh)

            for i, transcript in enumerate(db.children(gene, level=1), 1):

                if transcript.featuretype not in ["CDS", "exon"]:
                    update_transcript_attributes(gene, transcript, i)
                    print(transcript, file=fh)

                    for exon_or_cds in db.children(transcript, level=1):
                        print(exon_or_cds, file=fh)

                else:
                    # CDS/exon are a direct child of a gene; infer a transcript feature
                    inferred_transcript = infer_transcript(gene, i)
                    print(inferred_transcript, file=fh)

                    # now walk ALL children of the gene, update their Parent, and print
                    for exon_or_cds in db.children(gene, level=1):
                        exon_or_cds.attributes["Parent"] = inferred_transcript.attributes["ID"]
                        print(exon_or_cds, file=fh)
