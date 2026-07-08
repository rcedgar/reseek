#!/bin/bash -e

truth=sf

outdir=../big_scop40x
rm -rf $outdir/
mkdir -p $outdir
cd $outdir

db=../test_data/scop40x.bcb
lookup=../test_data/scop40x.lookup

for mode in fast sensitive
do
	name=$truth.$mode
	hits=$outdir/$name.hits
	reseek \
		-search $db \
		-$mode \
		-stats $truth \
		-columns query+target+pvalue \
		-db $db \
		-output $hits \
		-log $name.search.log

	reseek \
		-fast_bench_hits $hits \
		-lookup $lookup \
		-truth $truth \
		-log $name.sum3.log

	reseek \
		-fast_bench_hits $hits \
		-lookup $lookup \
		-truth top$truth \
		-log $name.top3.log
done
