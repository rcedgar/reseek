#!/usr/bin/python3

import sys

fn = sys.argv[1]

f = open(fn)

def get_flds(nf = None):
    line = f.readline()
    assert len(line) > 0
    flds = line[:-1].split('\t')
    if not nf is None:
        assert len(flds) == nf
    return flds

def assert_symmetrical(m):
    n = len(m[0])
    for i in range(n):
        assert len(m[i]) == n
        for j in range(i):
            d = m[i][j] - m[i][j]
            assert d == 0

def do_next():
    line = f.readline()
    if len(line) == 0:
        return None, None, None
    flds = line[:-1].split('\t')
    assert len(flds) == 4
    assert flds[0] == "FEATURE"
    name = flds[1]
    n = int(flds[2])
    ES = float(flds[3])

    S_ijs = []
    for i in range(n):
        row = []
        flds = get_flds(n+2)
        assert flds[0] == "S_ij"
        assert int(flds[1]) == i
        for j in range(n):
            S_ij = int(flds[j+2])
            row.append(S_ij)
        S_ijs.append(row)
    
    assert_symmetrical(S_ijs)

    return name, S_ijs, ES

names = []
ns = []
while True:
    name, S_ijs, ES = do_next()
    if name is None:
        break

    n = len(S_ijs[0])
    names.append(name)
    ns.append(n)

    print("\n//////////////////////////////////")
    print("// %s ES=%.4f" % (name, ES))
    print("//////////////////////////////////")

    print("int8_t %s_S_ij_vec[%d] = {" % (name, n*n))
    for i in range(n):
        s = ""
        for j in range(n):
            s += "%3s," % ("%d" % S_ijs[i][j])
        s += " // %2d" % i
        print(s)
    print(" };")

    print("int8_t %s_S_ij[%d][%d] = {" % (name, n, n))
    s =   "// 0 "
    for j in range(1, n):
        s += "%3d " % j
    print(s)
    for i in range(n):
        s = "{"
        for j in range(n):
            s += "%3s," % ("%d" % S_ijs[i][j])
        s += " }, // %2d" % i
        print(s)
    print("};")
