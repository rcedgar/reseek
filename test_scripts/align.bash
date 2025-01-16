#!/bin/bash -e

cd ../test_output

$reseek \
	-selfsearch dir.bca \
	-fast \
	-aln all_vs_all.aln \
	-log all_vs_all.log
