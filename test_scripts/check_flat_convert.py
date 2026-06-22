#!/usr/bin/env python3

import os
import sys

outdir = "../test_output/flat_convert/"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fasta import ReadSeqsDict

errors = 0
CCoordTol = 0.05

def err(msg):
    global errors
    errors += 1
    print("ERROR: %s" % msg, file=sys.stderr)

def read_cal_chains(fn):
    chains = {}
    label = None
    rows = []
    for line in open(fn):
        line = line.rstrip("\n\r")
        if not line:
            continue
        if line[0] == ">":
            if label is not None:
                chains[label] = rows
            label = line[1:]
            rows = []
        else:
            flds = line.split("\t")
            if len(flds) != 4:
                err("%s: expected 4 CAL fields, got %u in '%s'" % (fn, len(flds), line))
                continue
            aa = flds[0]
            if len(aa) != 1:
                err("%s: expected one aa letter, got '%s'" % (fn, aa))
                continue
            x, y, z = float(flds[1]), float(flds[2]), float(flds[3])
            rows.append((aa, x, y, z))
    if label is not None:
        chains[label] = rows
    return chains


def read_can_chains(fn):
    chains = {}
    label = None
    rows = []
    for line in open(fn):
        line = line.rstrip("\n\r")
        if not line:
            continue
        if line[0] == ">":
            if label is not None:
                chains[label] = rows
            label = line[1:]
            rows = []
        else:
            flds = line.split("\t")
            if len(flds) != 5:
                err("%s: expected 5 CAN fields, got %u in '%s'" % (fn, len(flds), line))
                continue
            aa = flds[0]
            if len(aa) != 1:
                err("%s: expected one aa letter, got '%s'" % (fn, aa))
                continue
            x, y, z = float(flds[1]), float(flds[2]), float(flds[3])
            hx = flds[4]
            if len(hx) != 2:
                err("%s: expected 2-digit hex nu, got '%s'" % (fn, hx))
                continue
            nu = int(hx, 16)
            if nu < 0 or nu > 255:
                err("%s: nu out of range: %s" % (fn, hx))
                continue
            rows.append((aa, x, y, z, nu))
    if label is not None:
        chains[label] = rows
    return chains


def decode_hexfasta(seq):
    if len(seq) % 2 != 0:
        return None
    out = []
    for i in range(0, len(seq), 2):
        out.append(int(seq[i:i + 2], 16))
    return out


def check_fasta_equal(fn_a, fn_b):
    a = ReadSeqsDict(fn_a)
    b = ReadSeqsDict(fn_b)
    if set(a.keys()) != set(b.keys()):
        err("%s and %s: label sets differ: %s vs %s" % (
            fn_a, fn_b, sorted(a.keys()), sorted(b.keys())))
        return
    for label in a:
        if a[label] != b[label]:
            err("%s and %s: sequence differs for %s" % (fn_a, fn_b, label))


def check_golden_fasta(out_fa, golden_fa):
    if not os.path.isfile(golden_fa):
        return
    check_fasta_equal(out_fa, golden_fa)


def check_cal_can_coords(cal_chains, can_chains):
    if set(cal_chains.keys()) != set(can_chains.keys()):
        err("CAL and CAN label sets differ")
        return
    for label in cal_chains:
        cal_rows = cal_chains[label]
        can_rows = can_chains.get(label)
        if can_rows is None:
            err("label %s missing from CAN" % label)
            continue
        if len(cal_rows) != len(can_rows):
            err("label %s: CAL len %u != CAN len %u" % (
                label, len(cal_rows), len(can_rows)))
            continue
        for i, (cal_row, can_row) in enumerate(zip(cal_rows, can_rows)):
            caa, cx, cy, cz = cal_row
            na, nx, ny, nz, _nu = can_row
            if caa != na:
                err("label %s pos %u: aa CAL=%s CAN=%s" % (label, i, caa, na))
            for v1, v2, axis in ((cx, nx, "x"), (cy, ny, "y"), (cz, nz, "z")):
                if abs(v1 - v2) > CCoordTol:
                    err("label %s pos %u: %s CAL=%.3f CAN=%.3f" % (
                        label, i, axis, v1, v2))


def check_nu_can_hex(can_chains, hex_fa):
    hexdict = ReadSeqsDict(hex_fa)
    if set(can_chains.keys()) != set(hexdict.keys()):
        err("CAN and hexfasta label sets differ")
        return
    for label in can_chains:
        nu_can = [row[4] for row in can_chains[label]]
        nu_hex = decode_hexfasta(hexdict[label])
        if nu_hex is None:
            err("label %s: odd-length hexfasta sequence" % label)
            continue
        if nu_can != nu_hex:
            err("label %s: CAN nu != hexfasta nu" % label)


def check_lengths(fa, hex_fa, kappa_fa):
    fadict = ReadSeqsDict(fa)
    hexdict = ReadSeqsDict(hex_fa)
    kappadict = ReadSeqsDict(kappa_fa)
    labels = set(fadict.keys())
    if labels != set(hexdict.keys()) or labels != set(kappadict.keys()):
        err("fasta / hexfasta / kappafasta label sets differ")
        return
    for label in labels:
        L = len(fadict[label])
        if len(hexdict[label]) != 2 * L:
            err("label %s: hexfasta len %u != 2*L=%u" % (
                label, len(hexdict[label]), 2 * L))
        if len(kappadict[label]) != L:
            err("label %s: kappafasta len %u != L=%u" % (
                label, len(kappadict[label]), L))
        if L == 0:
            err("label %s: zero-length sequence" % label)


def check_kappa_alphabet(kappa_fa):
    kappadict = ReadSeqsDict(kappa_fa)
    for label, seq in kappadict.items():
        for c in seq:
            if c < "A" or c > "z" or c in " \t":
                err("label %s: unexpected kappa char '%s'" % (label, c))


# outdir = sys.argv[1]
palms_fa = os.path.join(outdir, "palms.fa")
palms_cal = os.path.join(outdir, "palms.cal")
palms_can = os.path.join(outdir, "palms.can")
palms_hex = os.path.join(outdir, "palms.nu.hexfa")
palms_kappa = os.path.join(outdir, "palms.kappa.fa")
from_bcb = os.path.join(outdir, "from_bcb.fa")
from_can = os.path.join(outdir, "from_can.fa")
from_cal = os.path.join(outdir, "from_cal.fa")
golden = os.path.join(os.path.dirname(outdir), "..", "test_data", "mini.fa")
golden = os.path.normpath(golden)

for fn in (palms_fa, palms_cal, palms_can, palms_hex, palms_kappa,
            from_bcb, from_can, from_cal):
    if not os.path.isfile(fn):
        err("missing output file: %s" % fn)

if errors:
    print("%s FAILED: %u error(s)" %  (sys.argv[0], errors), file=sys.stderr)
    sys.exit(1)

check_fasta_equal(palms_fa, from_bcb)
check_fasta_equal(palms_fa, from_can)
check_fasta_equal(palms_fa, from_cal)
check_golden_fasta(palms_fa, golden)

cal_chains = read_cal_chains(palms_cal)
can_chains = read_can_chains(palms_can)
check_cal_can_coords(cal_chains, can_chains)
check_nu_can_hex(can_chains, palms_hex)
check_lengths(palms_fa, palms_hex, palms_kappa)
check_kappa_alphabet(palms_kappa)

if errors:
    print("%s FAILED: %u error(s)" %  (sys.argv[0], errors), file=sys.stderr)
    sys.exit(1)

print("%s SUCCESS" % sys.argv[0])
