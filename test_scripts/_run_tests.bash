#!/bin/bash -e

echo
echo _run_tests.bash
echo

START_SECS=$SECONDS

if [ x$reseek == x ] ; then
	echo "reseek not set"
	exit 1
fi

if [ ! -d ../big_data ] ; then
	echo "../big_data not found"
	exit 1
fi

if [ `uname -n` != "rip" ] ; then
	echo "Must run on rip (time and memory use checks depend on PC)"
	exit 1
fi

rm -rf ../test_output
mkdir ../test_output
mkdir -p ../test_results

log=../test_output/TEST_LOG.txt

date=`date "+%Y-%m-%d/%H:%M:%S"`
ver=`$reseek --version | tr -d ' \n\r'`
echo $date $ver STARTED >> $log
git status >> $log

echo STARTED `date` >> $log

./flat_convert.bash
./scop40x.bash

python3 ./check_scop40x.py >> $log
python3 ./check_flat_convert.py >> $log


python3 ./update_success_list.py $ver $date

END_SECS=$SECONDS
HOURS=$((ELAPSED_TOTAL / 3600))
MINUTES=$(((ELAPSED_TOTAL % 3600) / 60))
SECONDS_LEFT=$((ELAPSED_TOTAL % 60))

HH=$(printf "%02d" $HOURS)
MM=$(printf "%02d" $MINUTES)
SS=$(printf "%02d" $SECONDS_LEFT)

echo "COMPLETED $date elapsed $HH:$MM:$SS" \
	| tee -a $log

echo "$date $ver SUCCESS elapsed $HH:$MM:$SS" >> $log
