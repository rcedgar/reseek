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
IN=../test_data/mini.cal

rm -rf "$OUT"
mkdir -p "$OUT"

echo "flat_convert smoke test: $IN"

"$reseek" \
	-flat_convert "$IN" \
	-cal "$OUT/all.cal" \
	-can "$OUT/all.can" \
	-bca "$OUT/all.bca" \
	-bcb "$OUT/all.bcb" \
	-fasta "$OUT/all.fa" \
	-hexfasta "$OUT/all.hex.fa" \
	-kappafasta "$OUT/all.kappa.fa" \
	-log "$OUT/convert.log"

"$reseek" \
	-flat_convert "$OUT/all.bcb" \
	-fasta "$OUT/from_bcb.fa" \
	-log "$OUT/from_bcb.log"

"$reseek" \
	-flat_convert "$OUT/all.can" \
	-fasta "$OUT/from_can.fa" \
	-log "$OUT/from_can.log"

"$reseek" \
	-flat_convert "$OUT/all.cal" \
	-fasta "$OUT/from_cal.fa" \
	-log "$OUT/from_cal.log"

python3 ./check_flat_convert.py "$OUT"
