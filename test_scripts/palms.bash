#!/bin/bash -e

q=../test_data/palms.bca
db=../test_data/scop40.bca

if [ -z "$reseek" ] ; then
	reseek=../bin/reseek
fi

mkdir -p ../test_output
cd ../test_output

for mode in fast sensitive verysensitive ; do
	$reseek \
		-search $q \
		-db $db \
		-columns query+target+pvalue+evalue+qlo+qhi+ql+tlo+thi+tl+qcovpct+tcovpct \
		-output palms_scop40.$mode.hits \
		-$mode \
		-log palms_scop40_$mode.log
done

for x in palm*.hits ; do
	echo
	echo === $x ===
	sort -gk3 $x | head
done
