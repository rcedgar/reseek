#!/usr/bin/python3

import sys

errors = 0

dom2sf = {}
for line in open("../test_data/scop40_sf.tsv"):
    flds = line[:-1].split('\t')
    dom = flds[0]
    sf = flds[1]
    dom2sf[dom] = sf

def readhits(fn, TP, TP1, FP1, FP):
    global errors, tp, tp1, fp, fp1
    tp = 0
    tp1 = 0
    fp = 0
    fp1 = 0
    for line in open(fn):
        flds = line[:-1].split('\t')
        if len(flds) != 3:
            errors += 1
            print("ERROR %s not 3 flds" % fn)
            return
        try:
            E = float(flds[0])
            dom1 = flds[1]
            dom2 = flds[2]
            sf1 = dom2sf[dom1]
            sf2 = dom2sf[dom2]
        except:
            errors += 1
            print("ERROR %s exception" % fn)
            return
        if dom1 == dom2:
            errors += 1
            print("ERROR %s dom1==dom2" % fn)
            return
        if sf1 == sf2:
            tp += 1
        else:
            fp += 1

        if E <= 1:
            if sf1 == sf2:
                tp1 += 1
            else:
                fp1 += 1
    if tp < TP*0.95:
        print("ERROR TP too low: TP=%d TP1=%d FP1=%d FP=%d %s" % (tp, tp1, fp1, fp, fn))
        sys.exit(1)
    if tp1 < TP1*0.98:
        print("ERROR TP1 too low: TP=%d TP1=%d FP1=%d FP=%d %s" % (tp, tp1, fp1, fp, fn))
        sys.exit(1)
    if fp1 > FP1*1.02:
        print("ERROR FP1 too high: TP=%d TP1=%d FP1=%d FP=%d %s" % (tp, tp1, fp1, fp, fn))
        sys.exit(1)

    print("ok TP=%d TP1=%d FP1=%d FP=%d %s" % (tp, tp1, fp1, fp, fn))

# ok TP=244536 TP1=136178 FP1=21240 FP=3550160 ../test_output/scop40.tsv
# ok TP=114662 TP1=114662 FP1=9259 FP=9259 ../test_output/scop40-fast.tsv
# ok TP=133114 TP1=133114 FP1=19087 FP=19087 ../test_output/scop40-evalue1.tsv

readhits("../test_output/scop40.tsv", TP=244536, TP1=136178, FP1=21240, FP=3550160)
readhits("../test_output/scop40-fast.tsv", TP=114662, TP1=114662, FP1=9259, FP=9259)
readhits("../test_output/scop40-evalue1.tsv", TP=133114, TP1=133114, FP1=19087, FP=19087)
