#!/bin/bash -e

un=`uname -o`
if [ "$un" == Cygwin ] ; then
	os=win
elif [ "$un" == GNU/Linux ] ; then
	os=linux
else
	echo Bad un=$un
	exit 1
fi

outdir=../big_scop40x
#####################
rm -rf $outdir/*.$os.*
##############
mkdir -p $outdir
cd $outdir

db=../test_data/scop40x.bcb
lookup=../test_data/scop40x.lookup

for mode in kappa # all
do
	for truth in fold sf fam
	do
		name=$truth.$mode.$os
		hits=$outdir/$name.hits
		reseek \
			-search $db \
			-$mode \
			-stats $truth \
			-db $db \
			-output $hits \
			-log $name.search.log

		reseek \
			-fast_bench_hits $hits \
			-lookup $lookup \
			-truth $truth \
			-log $name.sum3.log

		if [ $truth != fam ] ; then
			reseek \
				-fast_bench_hits $hits \
				-lookup $lookup \
				-truth top$truth \
				-log $name.top3.log
		fi
	done
done

grep 3= $outdir/*3.log
