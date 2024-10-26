#!/bin/bash -e

cd ../test_output

$reseek \
	-search scop40.bca \
	-output scop40-fast.tsv \
	-fast \
	-log scop40-fast.log

$reseek \
	-search scop40.bca \
	-output scop40-sensitive.tsv \
	-sensitive \
	-log scop40-sensitive.log

$reseek \
	-search scop40.bca \
	-output scop40-evalue1.tsv \
	-fast \
	-evalue 1 \
	-log scop40-evalue1.log
