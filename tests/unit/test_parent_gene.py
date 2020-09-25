"""tests of parent gene assignment"""

import sys
import pytest
sys.path = ["..", "../.."] +  sys.path
from cat_test import get_input_file, get_output_file, diff_expected
from cat.parent_gene_assignment import assign_parents


@pytest.mark.parametrize("gene_set", ["HIGD1B", "DHX8", "TBCE", "PRR5-ARHGAP8"])
def test_parent(gene_set):
    """parent assignment cases from T2T CHM13 20200727 assembly"""

    parentsDf = assign_parents(get_input_file(f"{gene_set}.CHM13.tm-filtered.gp"),
                               get_input_file(f"{gene_set}.CHM13.tm-unfiltered.gp"),
                               get_input_file(f"CHM13.chrom.sizes"),
                               get_input_file(f"{gene_set}.CHM13.augPB.gp"))
    # group by gene id and then into in predictable order
    parentsDf.sort_values(["AssignedGeneId", "TranscriptId"], inplace=True)
    parentsDf.to_csv(get_output_file(".tsv"), sep='\t')
    diff_expected(".tsv")
