"""
This hacky script makes the LV3 FANTOM GTF work with CAT.
"""


def parse_gtf_attr_line(attr_line):
    """parse a GTF attributes line"""
    attr_line = [x.split(" ") for x in re.split("; +", attr_line.replace('"', ""))]
    attr_line[-1][-1] = attr_line[-1][-1].rstrip().replace(";", "")
    return dict(attr_line)


import re

recs = [x.rstrip().split("\t") for x in open("FANTOM_CAT.lv3_robust.gtf")]

attr_map = {}
new_recs = []
for rec in recs[1:]:
    rec = rec[:]
    attrs = parse_gtf_attr_line(rec[-1])
    if rec[2] == "gene":
        if "gene_id" not in attrs:
            attrs["gene_id"] = attrs["gene_name"]
        attrs["gene_biotype"] = attrs["geneSuperClass"]
        if attrs["gene_biotype"] == "all_mRNA":
            attrs["gene_biotype"] = "protein_coding"
        attr_map[attrs["gene_id"]] = attrs
        attrs["ID"] = attrs["gene_id"]
    elif rec[2] == "transcript":
        gene_attrs = attr_map[attrs["gene_id"]]
        attrs.update(gene_attrs)
        attrs["transcript_biotype"] = "protein_coding" if attrs["coding_status"] == "coding" else "non_coding"
        attrs["transcript_name"] = attrs["transcript_id"]
        attrs["Parent"] = attrs["gene_id"]
        attrs["ID"] = attrs["transcript_id"]
    elif rec[2] == "CDS" or rec[2] == "exon":
        attrs["Parent"] = attrs["transcript_id"]
    rec[-1] = ";".join(
        ["{}={}".format(x[0].lower() if x[0] != "Parent" and x[0] != "ID" else x[0], x[1]) for x in attrs.iteritems()]
    )
    new_recs.append(rec)


fh = open("FANTOM_CAT.lv3_robust.gff3", "w")
fh.write(recs[0][0] + "\n")
for rec in new_recs:
    fh.write("\t".join(rec) + "\n")


fh.close()
