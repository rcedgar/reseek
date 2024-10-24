#!/bin/bash -e

echo
echo _run_tests.bash
echo

cd ../src
rm -rf o/ ../bin/reseek*
./build_linux_x86.bash
cd ../test_scripts

rm -rf ../test_output
mkdir ../test_output
mkdir -p ../test_results

log=../test_output/TEST_LOG.txt

date=`date "+%Y-%m-%d/%H:%M:%S"`
ver=`../bin/reseek --version | tr -d ' \n\r'`
echo $date $ver STARTED >> $log
git status >> $log

echo STARTED `date` >> $log

./convert.bash
./align.bash
./columns.bash
./search.bash
./scop40.bash
./pdb2mega.bash

python3 ./check_logs.py >> $log
python3 ./check_convert.py >> $log
python3 ./check_columns.py >> $log
python3 ./check_scop40.py >> $log

python3 ./update_success_list.py $ver $date

echo COMPLETED $date >> $log

echo $date $ver SUCCESS >> $log
