#!/usr/bin/python3

errors = 0

dom2sf = {}
for line in open("../test_data/scop40_sf.tsv"):
    flds = line[:-1].split('\t')
    dom = flds[0]
    sf = flds[1]
    dom2sf[dom] = sf

def readhits(fn):
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
    print("ok TP=%d TP1=%d FP1=%d FP=%d %s" % (tp, tp1, fp1, fp, fn))
    return tp, tp1, fp

readhits("../test_output/scop40.tsv")
readhits("../test_output/scop40-fast.tsv")
readhits("../test_output/scop40-evalue1.tsv")
