#!/bin/bash -e

cd ../test_output

../bin/reseek \
	-search dir.bca \
	-aln all_vs_all.aln \
	-log all_vs_all.log
