#!/bin/bash -e

cd ../test_output

$reseek \
	-alignpair 5u0c.pdb \
	-input2 5wz3.pdb \
	-aln 5u0c_5wz3.pdb
	-log alignpair.log
