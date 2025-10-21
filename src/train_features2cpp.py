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
            d = abs(m[i][j] - m[i][j])
            assert d < 1e-6

def do_next():
    line = f.readline()
    if len(line) == 0:
        return None, None, None, None, None
    flds = line[:-1].split('\t')
    assert len(flds) == 4
    assert flds[0] == "FEATURE"
    name = flds[1]
    n = int(flds[2])
    ES = float(flds[3])

    f_is = []
    for i in range(n):
        flds = get_flds(3)
        assert flds[0] == "f_i"
        assert int(flds[1]) == i
        f_i = float(flds[2])
        f_is.append(f_i)

    f_ijs = []
    for i in range(n):
        row = []
        flds = get_flds(n+2)
        assert flds[0] == "f_ij"
        assert int(flds[1]) == i
        for j in range(n):
            f_ij = float(flds[j+2])
            row.append(f_ij)
        f_ijs.append(row)

    S_ijs = []
    for i in range(n):
        row = []
        flds = get_flds(n+2)
        assert flds[0] == "S_ij"
        assert int(flds[1]) == i
        for j in range(n):
            S_ij = float(flds[j+2])
            row.append(S_ij)
        S_ijs.append(row)
    
    assert_symmetrical(f_ijs)
    assert_symmetrical(S_ijs)

    return name, f_is, f_ijs, S_ijs, ES

print('#include "myutils.h"')
print('#include "features.h"')
print('#include "dss.h"')

names = []
ns = []
while True:
    name, f_is, f_ijs, S_ijs, ES = do_next()
    if name is None:
        break

    n = len(f_is)
    names.append(name)
    ns.append(n)

    print("\n//////////////////////////////////")
    print("// %s ES=%.4f" % (name, ES))
    print("//////////////////////////////////")
    print("\nstatic double %s_f_i[%d] = {" % (name, n))
    for i in range(n):
        print(f"   %.4f, // %d" % (f_is[i], i))
    print(" };")
    
    print("\nstatic double %s_f_ij[%d][%d] = {" % (name, n, n))
    for i in range(n):
        s = "   {"
        for j in range(n):
            s += " %10.4g," % f_ijs[i][j]
        s += " }, // %d" % i
        print(s)
    print(" };")

    print("\nstatic double %s_S_ij[%d][%d] = {" % (name, n, n))
    for i in range(n):
        s = "   {"
        for j in range(n):
            s += " %10.4g," % S_ijs[i][j]
        s += " }, // %d" % i
        print(s)
    print(" };")

print("")
print("float **g_ScoreMxs2[FEATURE_COUNT];")
print("float **g_FreqMxs2[FEATURE_COUNT];")
print("float *g_FreqVecs2[FEATURE_COUNT];")
print("uint g_AlphaSizes2[FEATURE_COUNT];");
print("")
print("static bool Init()")
print("\t{")
K = len(names)
for k in range(K):
    name = names[k]
    n = ns[k]
    print("\tasserta(DSSParams::GetAlphaSize(FEATURE_%s) == %d);" % (name, n))
print("")

for k in range(K):
    name = names[k]
    n = ns[k]
    print("\tg_AlphaSizes2[FEATURE_%s] = %d;" % (name, n))
print("")

for k in range(K):
    name = names[k]
    n = ns[k]
    print("\tg_FreqMxs2[FEATURE_%s] = myalloc(float *, %d);" % (name, n))
    print("\tg_ScoreMxs2[FEATURE_%s] = myalloc(float *, %d);" % (name, n))
    print("\tg_FreqVecs2[FEATURE_%s] = myalloc(float, %d);" % (name, n))
    print("\tfor (uint i = 0; i < %d; ++i)" % n)
    print("\t\t{")
    print("\t\tg_FreqVecs2[FEATURE_%s][i] = (float) %s_f_i[i];" % (name, name))
    print("\t\tg_FreqMxs2[FEATURE_%s][i] = myalloc(float, %d);" % (name, n))
    print("\t\tg_ScoreMxs2[FEATURE_%s][i] = myalloc(float, %d);" % (name, n))
    print("\t\tfor (uint j = 0; j < %d; ++j)" % n)
    print("\t\t\t{")
    print("\t\t\tg_FreqMxs2[FEATURE_%s][i][j] = (float) %s_f_ij[i][j];" % (name, name))
    print("\t\t\tg_ScoreMxs2[FEATURE_%s][i][j] = (float) %s_S_ij[i][j];" % (name, name))
    print("\t\t\t}")
    print("\t\t}")

print("\treturn true;")
print("\t}")
print("static bool InitDone = Init();")
