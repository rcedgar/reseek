#!/bin/bash -e

cd ../test_output

../bin/reseek \
	-pdb2mega ../test_structures/4v40.cif \
	-output ../test_output/4v40.mega \
	-log pdb2mega.log
