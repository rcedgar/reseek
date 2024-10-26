#!/bin/bash -e

cd ../test_output

$reseek \
	-search dir.bca \
	-fast \
	-output columns_default.tsv \
	-log columns_default.log

$reseek \
	-search dir.bca \
	-fast \
	-columns evalue+query+target \
	-output columns_same_as_default.tsv \
	-log columns_same_as_default.log

$reseek \
	-search dir.bca \
	-fast \
	-columns std \
	-output columns_std.tsv \
	-log columns_std.log

$reseek \
	-search dir.bca \
	-fast \
	-columns query+target+qlo+qhi+ql+tlo+thi+tl+cigar+qrow+trow \
	-output columns_local_rows.tsv \
	-log columns_local_rows.log

$reseek \
	-search dir.bca \
	-fast \
	-columns query+target+qlo+qhi+ql+tlo+thi+tl+cigar+qrowg+trowg \
	-output columns_global_rows.tsv \
	-log columns_global_rows.log
