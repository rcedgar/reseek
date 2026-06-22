#!/usr/bin/python3

import sys
import fasta

errors = 0

d1 = fasta.ReadSeqsDict("../test_data/dir.mu.fa")
d2 = fasta.ReadSeqsDict("../test_output/dir.mu.fa")

letters = 0
diffs = 0
def check():
    global errors
    global diffs
    labels1 = set(d1.keys())
    labels2 = set(d2.keys())
    if labels1 != labels2:
        errors += 1
        print("%s: ERROR nr labels differ", sys.argv[0])
        return
    for label in labels1:
        seq1 = d1[label]
        seq2 = d2[label]
        L1 = len(seq1)
        L2 = len(seq2)
        letters += L1
        if L1 != L2:
            errors += 1
            print("%s: ERROR lengths %d %d >%s", \
                sys.argv[0], L1, L2, label)
            return
        for i in range(L1):
            if seq1[i] != seq2[i]:
                diffs += 1
    if diffs/letters > 0.01:
        errors += 1
        print("%s: ERROR %d / %u diffs", \
            sys.argv[0], diffs, letters)

if errors == 0:            
    print("ok check_mu")
exit(1 if errors > 0 else 0)
