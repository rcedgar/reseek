#!/bin/bash -e

cd ../test_output

for idx in idxq idxt ; do
	/bin/time \
		-f "$idx elapsedsecs=%e memkb=%M pctcpu=%P systemsecs=%S usersecs=%U" \
		-o 1hhs_pdb90_$idx.time \
	$reseek \
		-search ../big_data/1hhs.pdb \
		-db ../big_data/pdb90.bca \
		-dbmu ../big_data/pdb90.mu.fa \
		-output 1hhs_pdb90_$idx.tsv \
		-fast \
		-$idx \
		-threads 1 \
		-log 1hhs_pdb90_$idx.log
done
