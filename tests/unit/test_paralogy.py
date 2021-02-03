""" tests for paralogous genes """ 

import sys
import pytest
sys.path = ["..", "../.."] +  sys.path
from cat_test import get_input_file, get_output_file, diff_expected
from cat.filter_transmap import filter_transmap
import luigi

@pytest.mark.parametrize("gene_set", ["TPS-all"])
# @pytest.mark.parametrize("gene_set", ["chr16"])
# @pytest.mark.parametrize("gene_set", ["whole-genome"])
def test_paralog(gene_set):
    """ paralogous gene test cases from T2T CHM13 full draft assembly """ 
    psl_target = luigi.LocalTarget(f"{gene_set}.CHM13.out.psl")
    json_target = luigi.LocalTarget(f"{gene_set}.CHM13.out.json")
    filteredDf = filter_transmap(get_input_file(f"{gene_set}.CHM13.psl"),
                                get_input_file("gencode.v35.annotation.gff3.psl"),
                                get_input_file(f"{gene_set}.CHM13.tm-unfiltered.gp"),
                                get_input_file("GRCh38.db"),
                                psl_tgt=psl_target,
                                global_near_best=0.15,    #default
                                # global_near_best=0.05,
                                filter_overlapping_genes=False,
                                # overlapping_ignore_bases=0,
                                overlapping_gene_distance=0,
                                json_tgt=json_target,
                                annotation_gp=get_input_file("gencode.v35.annotation.gff3.gp")) #json_target
    filteredDf.to_csv(get_output_file(".tsv"), sep='\t')
    # filteredDf.sort_values(["AssignedGeneId", "TranscriptId"], inplace=True)
    diff_expected(".tsv")

