#!/bin/bash -e

reseek=../github_releases/reseek-v2.8-linux-x86

if [ ! -x $reseek ] ; then
	echo Download ../github_releases/reseek-v2.8-linux-x86 and chmod +x
	exit 1
fi

truth=sf

outdir=../big_scop40x_v2
rm -rf $outdir/
mkdir -p $outdir
cd $outdir

db=../test_data/scop40x.bca
lookup=../test_data/scop40x.lookup

for mode in fast sensitive
do
	name=v2_$mode
	hits=$outdir/$name.hits
	$reseek \
		-search $db \
		-$mode \
		-columns query+target+evalue \
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
