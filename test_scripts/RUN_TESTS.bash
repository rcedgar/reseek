#!/bin/bash -e

cd ../src
./build_linux.py
cd ../test_scripts

rm -rf ../test_output
mkdir ../test_output
mkdir -p ../test_results

date=`date "+%Y-%m-%d/%H:%M:%S"`
ver=`reseek --version`
echo $date $ver INCOMPLETE/FAILED \
  > ../test_results/test_result.txt

log=../test_output/TEST_LOG.txt

echo STARTED `date` | tee $log

./convert.bash
./align.bash
./columns.bash
./search.bash
./scop40.bash

./check_logs.py | tee -a $log
./check_convert.py | tee -a $log
./check_columns.py | tee -a $log
./check_scop40.py | tee -a $log

./update_success_list.py $ver $date

echo COMPLETED `date` | tee -a $log

echo $date "$ver" SUCCESS \
  > ../test_results/test_result.txt
