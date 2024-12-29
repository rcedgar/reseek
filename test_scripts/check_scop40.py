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
        if len(flds) != 4:
            errors += 1
            print("%s: ERROR %s not 3 flds" % (sys.argv[0], fn))
            return
        try:
            E = float(flds[3])
            dom1 = flds[1]
            dom2 = flds[2]
            sf1 = dom2sf[dom1]
            sf2 = dom2sf[dom2]
        except:
            errors += 1
            print("%s: ERROR %s exception" % (sys.argv[0], fn))
            return
        if dom1 == dom2:
            continue
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
        print("ERROR TP too low: TP=%d, TP1=%d, FP1=%d, FP=%d %s" % (tp, tp1, fp1, fp, fn))
        sys.exit(1)
    if tp1 < TP1*0.98:
        print("ERROR TP1 too low: TP=%d, TP1=%d, FP1=%d, FP=%d %s" % (tp, tp1, fp1, fp, fn))
        sys.exit(1)
    if fp1 > FP1*1.02:
        print("ERROR FP1 too high: TP=%d, TP1=%d, FP1=%d, FP=%d %s" % (tp, tp1, fp1, fp, fn))
        sys.exit(1)
    if fp > FP*1.02:
        print("ERROR FP too high: TP=%d, TP1=%d, FP1=%d, FP=%d %s" % (tp, tp1, fp1, fp, fn))
        sys.exit(1)

    print("ok TP=%d, TP1=%d, FP1=%d, FP=%d %s" % (tp, tp1, fp1, fp, fn))

readhits("../test_output/scop40-fast.tsv", TP=175928, TP1=133730, FP1=8300, FP=125400)
readhits("../test_output/scop40-sensitive.tsv", TP=223980, TP1=148640, FP1=8750, FP=161000)
readhits("../test_output/scop40-evalue1.tsv", TP=132754, TP1=132754, FP1=7800, FP=7800)
