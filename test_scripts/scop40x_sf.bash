#!/bin/bash -e

START_SECS=$SECONDS

un=`uname -o`
if [ "$un" == Cygwin ] ; then
	os=win
elif [ "$un" == GNU/Linux ] ; then
	os=linux
else
	echo Bad un=$un
	exit 1
fi

truth=sf

outdir=../big_scop40x
#####################
rm -rf $outdir/$truth.*.$os.*
##############
mkdir -p $outdir
cd $outdir

db=../test_data/scop40x.bcb
lookup=../test_data/scop40x.lookup

for mode in all kappa
do
	name=$truth.$mode.$os
	hits=$outdir/$name.hits
	reseek \
		-flat_search_$mode $db \
		-varstr =$truth \
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

python3 ../test_scripts/check_scop40x.py sf

END_SECS=$SECONDS
ELAPSED_TOTAL=$((END_SECS - START_SECS))
HOURS=$((ELAPSED_TOTAL / 3600))
MINUTES=$(((ELAPSED_TOTAL % 3600) / 60))
SECONDS_LEFT=$((ELAPSED_TOTAL % 60))

HH=$(printf "%02d" $HOURS)
MM=$(printf "%02d" $MINUTES)
SS=$(printf "%02d" $SECONDS_LEFT)

echo "elapsed $HH:$MM:$SS"
