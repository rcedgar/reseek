#!/bin/bash -e

cd ../test_output

$reseek \
	-search dir.bca \
	-fast \
	-output search_dir.bca.tsv \
	-log search_dir_bca.log

$reseek \
	-search dir.cal \
	-output search_dir.cal.tsv \
	-fast \
	-log search_dir_cal.log

$reseek \
	-search dir.cal \
	-db dir.bca \
	-output search_dir_cal_bca.tsv \
	-fast \
	-log search_dir_cal_bca.log

$reseek \
	-search dir.bca \
	-db dir.cal \
	-output search_dir_bca_cal.tsv \
	-fast \
	-log search_dir_bca_cal.log
