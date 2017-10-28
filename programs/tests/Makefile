##
# - To run ncbi_submit_tests and pass through table2asn, you must first get
#   the fasta file for the test chromosome sequence and place in this directory
#   under the name chimp1.fa.  This file is too big to check in.
#   FIXME: use git large file support
#

cat_to_ncbi_submit = ../cat_to_ncbi_submit
tblToDumbBed = bin/tblToDumbBed

ncbibin = ../../../ncbi/$(shell uname -s)/bin
table2asn = $(ncbibin)/table2asn
tbl2asn =  $(ncbibin)/tbl2asn

# control running using tbl2asn, table2asn, or both
# tbl2asn is an option because the current version of
# table2asn will not generate genbank flat files
useTable2asn = no
useTbl2asn = yes

# overide setting based on the existance of the programs and the
# test fasta
ifeq ($(wildcard chimp1.fa),)
  useTable2asn = no
  useTbl2asn = no
endif
ifeq ($(wildcard ${table2asn}),)
  useTable2asn = no
endif
ifeq ($(wildcard ${tbl2asn}),)
  useTbl2asn = no
endif

diff = diff -u

test: ncbi_submit_tests


ncbi_submit_tests: ncbi_submit_test1



##
# tests on V1 pacbio
# /hive/groups/recon/projs/primates/susie_the_orangutan/annotation_fix_indels/
#      genome_files/Clint_Chimp.fa
#      consensus_gene_set/Clint_Chimp.fixed.gp_info
#      consensus_gene_set/Clint_Chimp.gp
#
# Tests are on sequence 000000F_1_80365438_quiver_pilon, which should be in a file chimp1.fa in tests directory.
#
# Ones with in-frame stops even after frame-indel adjustment:
#    Clint_Chim_T0000047 - pos strand, deletion in first coding exon, should drop as too short
#    Clint_Chim_T0000217 - pos strand, deletion near stop
#    Clint_Chim_T0000883 - neg strand, single exon near start
#    Clint_Chim_T0001227 - neg strand, near stop
# 

ncbi_submit_test1: mkdirs
	${cat_to_ncbi_submit} input/chimp1.gp input/chimp1.gp_info fred output/$@.tbl
	${diff} expected/$@.tbl output/$@.tbl
	${tblToDumbBed} output/$@.tbl output/$@.bed
ifeq (${useTable2asn}, yes)
	rm -rf output/$@.table2asn
	mkdir -p  output/$@.table2asn
	${table2asn} -V vbag -M n -j "[organism=Pan troglodytes]" -J -c - -euk -t input/chimp1.sbt -i chimp1.fa -f output/$@.tbl -o output/$@.table2asn/chimp1.sqn -Z output/$@.table2asn/chimp1.report
endif
ifeq (${useTbl2asn}, yes)
	rm -rf output/$@.tbl2asn
	mkdir -p  output/$@.tbl2asn
	${tbl2asn} -V vb -M n -c s -j "[organism=Pan troglodytes]" -t input/chimp1.sbt -i chimp1.fa -f output/$@.tbl -o output/$@.tbl2asn/chimp1.sqn -Z output/$@.tbl2asn/chimp1.report
	mv chimp1.gbf chimp1.val errorsummary.val output/$@.tbl2asn/
endif


mkdirs:
	@mkdir -p output

clean:
	rm -rf output
