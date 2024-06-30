#!/bin/bash -e

cd ../test_output

../bin/reseek \
	-search dir.bca \
	-output columns_default.tsv \
	-log columns_default.log

../bin/reseek \
	-search dir.bca \
	-columns evalue+query+target \
	-output columns_same_as_default.tsv \
	-log columns_same_as_default.log

../bin/reseek \
	-search dir.bca \
	-columns std \
	-output columns_std.tsv \
	-log columns_std.log

../bin/reseek \
	-search dir.bca \
	-columns query+target+qlo+qhi+ql+tlo+thi+tl+cigar+qrow+trow \
	-output columns_local_rows.tsv \
	-log columns_local_rows.log

../bin/reseek \
	-search dir.bca \
	-columns query+target+qlo+qhi+ql+tlo+thi+tl+cigar+qrowg+trowg \
	-output columns_global_rows.tsv \
	-log columns_global_rows.log
