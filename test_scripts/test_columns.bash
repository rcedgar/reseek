#!/bin/bash -e

truth=superfamily
mode=fast
db=../test_data/palms.bcb
fa=palms.fa
hits=../test_output/palms.hits
cols=query+target+qlo+qhi+ql+tlo+thi+tl+cigar+qrow+trow+pvalue

mkdir -p ../test_output

reseek \
	-flat_convert $db \
	-fasta $fa \
	-log test_columns_convert_fasta.log

reseek \
	-search $db \
	-$mode \
	-db $db \
	-stats $truth \
	-columns $cols \
	-output $hits \
	-log test_columns_search.log

#    fasta_fn = sys.argv[1]
#    columns = sys.argv[2].split("+")
#    tsv_fn = sys.argv[3]
python ../test_scripts/check_hits_against_fasta.py \
	$fa \
	$cols \
	$hits \
	| tee check_colums.result
