#!/bin/bash -e

cd ../test_output

$reseek \
	-alignpair ../test_data/5u0c.pdb \
	-input2 ../test_data/5wz3.pdb \
	-aln 5u0c_5wz3.aln \
	-output 5u0c_5wz3.pdb \
	-log alignpair.log
