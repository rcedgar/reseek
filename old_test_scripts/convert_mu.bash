#!/bin/bash -e

cd ../test_output

$reseek \
	-convert2mu dir.bca \
	-fasta dir.mu.fa \
	-log convert2mu.log
