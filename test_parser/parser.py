# test for gff parser
from pandas.util.testing import assert_frame_equal
import pandas as pd
import sys
import tools.gff3


def test_ncbi():
    test = tools.gff3.parse_gff3("test_parser/ncbiGp", "test_parser/ncbiAttrs")
    ncbi = pd.read_csv("test_parser/ncbiTest", sep="\t", index_col=0)
    assert_frame_equal(test, ncbi, check_dtype=False)


def test_gencode():
    test = tools.gff3.parse_gff3("test_parser/gencodeGp", "test_parser/gencodeAttrs")
    gencode = pd.read_csv("test_parser/gencodeTest", sep="\t", index_col=0)
    assert_frame_equal(test, gencode, check_dtype=False)


def test_ensembl():
    test = tools.gff3.parse_gff3("test_parser/ensemblGp", "test_parser/ensemblAttrs")
    ensembl = pd.read_csv("test_parser/ensemblTest", sep="\t", index_col=0)
    assert_frame_equal(test, ensembl, check_dtype=False)
