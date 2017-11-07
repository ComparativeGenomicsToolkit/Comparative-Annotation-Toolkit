##
# - To run ncbi_submit_tests and pass through table2asn, you must first get
#   the fasta file for the test chromosome sequence and place in this directory
#   under the name data/chimp1.fa.  These files are too big to check in.
#   They are generate on demand.
#
#   - ncbi_submit_test1 - use V2 of chimp.
#
#   FIXME: use git large file support for fasta files
#
SHELL = /bin/bash

cat_to_ncbi_submit = ../cat_to_ncbi_submit
tblToDumbBed = bin/tblToDumbBed

ncbibin = /hive/groups/recon/local/bin
tbl2asn =  $(ncbibin)/tbl2asn

chimpSrcTwoBit = /hive/groups/recon/projs/primates/primates_indel_corrected_bionano_cut/genome_files/Clint_Chimp.2bit

# tbl2asn is an option because the current version of
useTbl2asn = yes

# overide setting based on the existance of the programs and the
# test fasta
ifeq ($(wildcard data/chimp1.fa),)
  useTbl2asn = no
endif
ifeq ($(wildcard ${tbl2asn}),)
  useTbl2asn = no
endif

diff = diff -u

test: ncbi_submit_tests


ncbi_submit_tests: ncbi_submit_test1


ncbi_submit_test1: mkdirs data/chimp1.fa
	${cat_to_ncbi_submit} input/chimp1.gp input/chimp1.gp_info fred output/$@.tbl
	${diff} expected/$@.tbl output/$@.tbl
ifeq (${useTbl2asn}, yes)
	rm -rf output/$@.tbl2asn
	mkdir -p output/$@.tbl2asn
	ln -f data/chimp1.fa output/$@.tbl2asn/
	${tbl2asn} -V vb -M n -c s -j "[organism=Pan troglodytes]" -t input/chimp1.sbt -i output/$@.tbl2asn/chimp1.fa -f output/$@.tbl -o output/$@.tbl2asn/chimp1.sqn -Z output/$@.tbl2asn/chimp1.report
endif

# extract genome sequence from 2bit
data/chimp1.fa: ${chimpSrcTwoBit} input/chimp1.gp
	@mkdir -p data
	twoBitToFa -seqList=<(cut -f 2 input/chimp1.gp | sort -u) ${chimpSrcTwoBit} $@.tmp
	mv -f $@.tmp $@

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
