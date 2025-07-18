#!/bin/bash -e

cd ../test_output

$reseek \
	-search scop40.bca \
	-db scop40.bca \
	-output scop40-fast.tsv \
	-fast \
	-log scop40-fast.log

$reseek \
	-search scop40.bca \
	-db scop40.bca \
	-output scop40-sensitive.tsv \
	-sensitive \
	-log scop40-sensitive.log

if [ "$dovs" == no ] ; then
	echo =================================
	echo SKIPPING --verysensitive
	echo =================================
else
	$reseek \
		-search scop40.bca \
		-db scop40.bca \
		-output scop40-verysensitive.tsv \
		-verysensitive \
		-log scop40-verysensitive.log
fi

$reseek \
	-search scop40.bca \
	-db scop40.bca \
	-output scop40-evalue1.tsv \
	-fast \
	-evalue 1 \
	-log scop40-evalue1.log
