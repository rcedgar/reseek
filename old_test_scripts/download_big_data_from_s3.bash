#!/bin/bash -e

mkdir -p ../big_data
cd ../big_data

aws s3 cp  s3://serratus-rce-mirror/reseek/test_data/big_data.tz .

tar -zxvf big_data.tz
