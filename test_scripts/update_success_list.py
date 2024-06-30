#!/usr/bin/python3

import sys

ver = sys.argv[1]
date = sys.argv[2]

fn = "../test_results/success_list.txt"
lines = []
try:
    for line in open(fn):
        lines.append(line[:-1])
except:
    pass

found = False
for line in lines:
    flds = line.split('\t')
    if flds[0] == ver:
        print("Build previously suceeded %s" % line)
        sys.exit(0)

f = open(fn, "a")
line = ver + "\t" + date
f.write(line + "\n")
f.close()
print("Build suceeded %s" % line)
