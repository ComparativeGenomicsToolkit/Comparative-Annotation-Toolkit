##
# - To run ncbi_submit_tests and pass through table2asn, you must first get
#   the fasta file for the test chromosome sequence and place in this directory
#   under the name data/chimp{1,2}.fa.  These files are too big to check in.
#
#   - ncbi_submit_test1 - use V1 of chimp, which is not indel corrected.  This induces more
#     corner cases.  Required sequences in are in input/chimp1.chroms
#   - ncbi_submit_test2 - use V1 of chimp, which is not indel corrected.  This tests newer issues
#     Required sequences in are in input/chimp2.chroms
#
#   FIXME: use git large file support for fasta files
#

cat_to_ncbi_submit = ../cat_to_ncbi_submit
tblToDumbBed = bin/tblToDumbBed

ncbibin = /hive/groups/recon/local/bin
tbl2asn =  $(ncbibin)/tbl2asn

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


ncbi_submit_tests: ncbi_submit_test1 ncbi_submit_test2


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
ifeq (${useTbl2asn}, yes)
	rm -rf output/$@.tbl2asn
	mkdir -p output/$@.tbl2asn
	ln -f data/chimp1.fa output/$@.tbl2asn/
	${tbl2asn} -V vb -M n -c s -j "[organism=Pan troglodytes]" -t input/chimp1.sbt -i output/$@.tbl2asn/chimp1.fa -f output/$@.tbl -o output/$@.tbl2asn/chimp1.sqn -Z output/$@.tbl2asn/chimp1.report
endif


ncbi_submit_test2: mkdirs
	${cat_to_ncbi_submit} input/chimp2.gp input/chimp2.gp_info fred output/$@.tbl
	${diff} expected/$@.tbl output/$@.tbl
ifeq (${useTbl2asn}, yes)
	rm -rf output/$@.tbl2asn
	mkdir -p output/$@.tbl2asn
	ln -f data/chimp2.fa output/$@.tbl2asn/
	${tbl2asn} -V vb -M n -c s -j "[organism=Pan troglodytes]" -t input/chimp1.sbt -i output/$@.tbl2asn/chimp2.fa -f output/$@.tbl -o output/$@.tbl2asn/chimp2.sqn -Z output/$@.tbl2asn/chimp2.report
endif

updateFa: data/chimp2.fa

# FIXME: make this more stable; this is currently just the subset
chimpSrcFa = /hive/groups/recon/projs/primates/primates_indel_corrected_bionano_cut/genome_files/Clint_Chimp.fa

data/chimp2.fa: input/chimp2.chroms ${chimpSrcFa}
	faSomeRecords ${chimpSrcFa} input/chimp2.chroms $@.tmp
	mv -f $@.tmp $@

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
