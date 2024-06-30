#!/bin/bash -e

cd ../test_output

../bin/reseek \
	-search scop40.bca \
	-output scop40.tsv \
	-minchainlength 5 \
	-log scop40.log

../bin/reseek \
	-search scop40.bca \
	-output scop40-fast.tsv \
	-minchainlength 5 \
	-fast \
	-log scop40-fast.log

../bin/reseek \
	-search scop40.bca \
	-output scop40-evalue1.tsv \
	-minchainlength 5 \
	-evalue 1 \
	-log scop40-evalue1.log
