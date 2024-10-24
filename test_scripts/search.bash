#!/bin/bash -e

cd ../test_output

../bin/reseek \
	-search dir.bca \
	-output search_dir.bca.tsv \
	-log search_dir_bca.log

../bin/reseek \
	-search dir.cal \
	-output search_dir.cal.tsv \
	-log search_dir_cal.log

../bin/reseek \
	-search dir.cal \
	-db dir.bca \
	-output search_dir_cal_bca.tsv \
	-log search_dir_cal_bca.log

../bin/reseek \
	-search dir.bca \
	-db dir.cal \
	-output search_dir_bca_cal.tsv \
	-log search_dir_bca_cal.log
