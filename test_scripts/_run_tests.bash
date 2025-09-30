#!/bin/bash -e

echo
echo _run_tests.bash
echo

if [ x$reseek == x ] ; then
	echo "reseek not set"
	exit 1
fi

if [ ! -d ../big_data ] ; then
	echo "../big_data not found"
	exit 1
fi

if [ `uname -n` != "rip" ] ; then
	echo Must run on rip
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

./idxqt_speed.bash
./convert.bash
./align.bash
./alignpair.bash
./columns.bash
./search.bash
./scop40.bash
./pdb2mega.bash

python3 ./check_idxqt_speed.py >> $log
python3 ./check_logs.py >> $log
python3 ./check_convert.py >> $log
python3 ./check_columns.py >> $log
python3 ./check_scop40.py >> $log

python3 ./update_success_list.py $ver $date

echo COMPLETED $date >> $log

echo $date $ver SUCCESS >> $log
