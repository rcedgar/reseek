#!/bin/bash -e

cd ../src
./build_linux.py
cd ../test_scripts

rm -rf ../test_output
mkdir ../test_output
mkdir -p ../test_results

log=../test_output/TEST_LOG.txt

date=`date "+%Y-%m-%d/%H:%M:%S"`
ver=`reseek --version | tr -d ' \n\r'`
echo $date $ver STARTED >> $log
git status >> $log

echo STARTED `date` >> $log

./convert.bash
./align.bash
./columns.bash
./search.bash
./scop40.bash

./check_logs.py >> $log
./check_convert.py >> $log
./check_columns.py >> $log
./check_scop40.py >> $log

./update_success_list.py $ver $date

echo COMPLETED $date >> $log

echo $date $ver SUCCESS >> $log
