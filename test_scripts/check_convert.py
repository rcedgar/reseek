#!/usr/bin/python3

import sys

def read_fasta(fn):
    seqdict = {}
    label = None
    seq = ""
    for line in open(fn):
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == ">":
            if len(seq) > 0:
                seqdict[label] = seq
            label = line[1:]
            seq = ""
        else:
            seq += line
    seqdict[label] = seq
    return seqdict

errors = 0
def check_fasta(fn):
    global errors
    seqdict2 = read_fasta(fn)
    for label in chains.keys():
        seq = chains[label]
        try:
            seq2 = seqdict2[label]
        except:
            errors += 1
            print(seqdict2.keys())
            print("%s: ERROR %s not found in %s" % (sys.argv[0], label, fn))
            return
        if seq == seq2:
            if 0:
                print("seq ok %s" % label)
        else:
            errors += 1
            print("%s: ERROR %s different %s" % (sys.argv[0], label, fn))
            return
    print("ok %s" % fn)

chains = read_fasta("../test_data/chains.fa")
check_fasta("../test_output/dir.fa")
check_fasta("../test_output/files.fa")
check_fasta("../test_output/cvt_dir_cal.fa")

exit(1 if errors > 0 else 0)
