#!/bin/bash

if [ ! -s ./_run_tests.bash ] ; then
	echo "ERROR -- must run in test_scripts/ directory"
	exit 1
fi

if [ x$1 == x ] ; then
	./recompile.bash
	export reseek=../bin/reseek
else
	if [ ! -x $1 ] ; then
		echo Not executable $x
		exit 1
	fi
	export reseek=$1
fi

./_run_tests.bash

if [ $? != 0 ] ; then
	echo === FAILED ===
fi

tail ../test_output/TEST_LOG.txt
