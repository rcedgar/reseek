#!/usr/bin/python3

import os
import sys
import scop40

errors = 0

sc = scop40.Scop40('e', "sf2", "../test_data/dom_scopid.tsv")

def dofile(tsv_fn, SEPQ0_1, SEPQ1, SEPQ10):
	global errors

	name = tsv_fn.split('/')[-1].split('.')[0]

	qfldnr = 0
	tfldnr = 1
	scorefldnr = 2
	is_sorted = False
	sc.eval_file(tsv_fn, qfldnr, tfldnr, scorefldnr, is_sorted)
	sepq0_1 = sc.tpr_at_fpepq0_1
	sepq1 = sc.tpr_at_fpepq1
	sepq10 = sc.tpr_at_fpepq10

	d0_1 = sepq0_1-SEPQ0_1
	d1 = sepq1-SEPQ1
	d10 = sepq10-SEPQ10

	if d0_1 < -0.01:
		errors += 1
		print("%s: ERROR d0_1" % sys.argv[0])
	if d1 < -0.01:
		errors += 1
		print("%s: ERROR d1" % sys.argv[0])
	if d10 < -0.01:
		errors += 1
		print("%s: ERROR d10" % sys.argv[0])

	s = "SEPQ0.1=%.4f(%+.4f)" % (sepq0_1, d0_1)
	s += " SEPQ1=%.4f(%+.4f)" % (sepq1, d1)
	s += " SEPQ10=%.4f(%+.4f)" % (sepq10, d10)
	s += " " + name
	print(s)

# SEPQ0.1=0.2107 SEPQ1=0.3144 SEPQ10=0.3882 S1FP=0.3350 N1FP=152340 area=7.14 fast
# SEPQ0.1=0.2173 SEPQ1=0.3411 SEPQ10=0.4745 S1FP=0.3874 N1FP=176191 area=10.6 sensitive
# SEPQ0.1=0.2106 SEPQ1=0.2950 SEPQ10=0.2950 S1FP=0.2848 N1FP=129529 area=3.99 evalue1

dofile("../test_output/scop40-evalue1.tsv", 		0.2100, 0.3100, 0.3500)
dofile("../test_output/scop40-fast.tsv", 			0.2100, 0.3140, 0.4200)
dofile("../test_output/scop40-sensitive.tsv",		0.2170, 0.3410, 0.4740)

if errors == 0:
	print("%s: PASSED" % sys.argv[0])
exit(1 if errors > 0 else 0)
