#!/bin/bash

mkdir -p ../tmp

log=../tmp/TEST_LOG_NO_E.txt

python3 ./check_idxqt_speed.py >> $log
python3 ./check_logs.py >> $log
python3 ./check_convert.py >> $log
python3 ./check_columns.py >> $log
python3 ./check_scop40.py >> $log

echo ERROR `grep -c ERROR $log`
echo OK `grep -c ^ok $log`
ls -lh $log
