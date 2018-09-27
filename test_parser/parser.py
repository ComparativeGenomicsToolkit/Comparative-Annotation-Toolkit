# test for gff parser
from pandas.util.testing import assert_frame_equal
import pandas as pd

def test_ncbi():
   test = run("ncbiGp","ncbiAttrs")
   ncbi = pd.read_csv("ncbiTest",sep='\t',index_col=0)
   assert_frame_equal(test,ncbi,check_dtype=False)

def test_gencode():
   test = run("gencodeGp","gencodeAttrs")
   gencode = pd.read_csv("gencodeTest",sep='\t',index_col=0)
   assert_frame_equal(test,gencode,check_dtype=False)

def test_ensembl():
   test = run("ensemblGp","ensemblAttrs")
   ensembl = pd.read_csv("ensemblTest",sep='\t',index_col=0)
   assert_frame_equal(test,ensembl,check_dtype=False)

