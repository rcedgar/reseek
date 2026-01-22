#!/bin/bash -e

cd ../test_output

if [ -z "$reseek" ] ; then
	reseek=../bin/reseek
fi

mkdir -p ../test_output
cd ../test_output

$reseek \
	-search scop40.bca \
	-db scop40.bca \
	-output scop40-fast.tsv \
	-columns query+target+evalue \
	-fast \
	-log scop40-fast.log

$reseek \
	-search scop40.bca \
	-db scop40.bca \
	-output scop40-sensitive.tsv \
	-columns query+target+evalue \
	-sensitive \
	-log scop40-sensitive.log

$reseek \
	-search scop40.bca \
	-db scop40.bca \
	-output scop40-evalue1.tsv \
	-columns query+target+evalue \
	-fast \
	-evalue 1 \
	-log scop40-evalue1.log
