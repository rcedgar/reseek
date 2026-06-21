#!/bin/bash -e

if [ x"$reseek" == x ] ; then
	reseek=../bin/reseek
fi

if [ ! -x "$reseek" ] ; then
	echo "ERROR: reseek not executable: $reseek"
	echo "Set reseek=PATH or build ../bin/reseek"
	exit 1
fi

OUT=../test_output/flat_convert
IN=../test_data/palms.bca

rm -rf "$OUT"
mkdir -p "$OUT"

"$reseek" \
	-flat_convert "$IN" \
	-cal "$OUT/palms.cal" \
	-can "$OUT/palms.can" \
	-bca "$OUT/palms.bca" \
	-bcb "$OUT/palms.bcb" \
	-fasta "$OUT/palms.fa" \
	-nuhexfasta "$OUT/palms.nu.hexfa" \
	-kappafasta "$OUT/palms.kappa.fa" \
	-log "$OUT/convert.log"

"$reseek" \
	-flat_convert "$OUT/palms.bcb" \
	-fasta "$OUT/from_bcb.fa" \
	-log "$OUT/from_bcb.log"

"$reseek" \
	-flat_convert "$OUT/palms.can" \
	-fasta "$OUT/from_can.fa" \
	-log "$OUT/from_can.log"

"$reseek" \
	-flat_convert "$OUT/palms.cal" \
	-fasta "$OUT/from_cal.fa" \
	-log "$OUT/from_cal.log"

python3 ./check_flat_convert.py "$OUT"
