#!/bin/bash -e

cd ../test_output

../bin/reseek \
	-search scop40.bca \
	-output scop40.tsv \
	-log scop40.log

../bin/reseek \
	-search scop40.bca \
	-output scop40-sensitive.tsv \
	-sensitive \
	-log scop40-sensitive.log

../bin/reseek \
	-search scop40.bca \
	-output scop40-evalue1.tsv \
	-evalue 1 \
	-log scop40-evalue1.log
