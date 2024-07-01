#!/bin/bash -e

./_run_tests.bash

if [ $? != 0 ] ; then
	echo === FAILED ===
fi

tail ../test_output/TEST_LOG.txt
