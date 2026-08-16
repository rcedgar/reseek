#!/usr/bin/env python3
"""Synthetic CAL -> BCB with 1 vs N threads; round-trip coords must match."""

import os
import sys
import subprocess
import tempfile

CCoordTol = 0.05


def err(msg):
    print("ERROR: %s" % msg)
    raise SystemExit(1)


def write_cal(fn, chains):
    with open(fn, "w", newline="\n") as f:
        for label, rows in chains:
            f.write(">%s\n" % label)
            for aa, x, y, z in rows:
                f.write("%c\t%.1f\t%.1f\t%.1f\n" % (aa, x, y, z))


def make_chain(label, L, x0):
    aas = "ACDEFGHIKLMNPQRSTVWY"
    rows = []
    for i in range(L):
        aa = aas[i % len(aas)]
        rows.append((aa, x0 + i, float(i), float(-i)))
    return (label, rows)


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
            if len(flds) != 4 or len(flds[0]) != 1:
                err("%s: bad CAL line '%s'" % (fn, line))
            aa = flds[0]
            x, y, z = float(flds[1]), float(flds[2]), float(flds[3])
            rows.append((aa, x, y, z))
    if label is not None:
        chains[label] = rows
    return chains


def compare_cal(a, b, name_a, name_b):
    if set(a) != set(b):
        err("label set mismatch %s vs %s: %s vs %s" %
            (name_a, name_b, sorted(a), sorted(b)))
    for lab in a:
        ra, rb = a[lab], b[lab]
        if len(ra) != len(rb):
            err("%s L=%u vs %s L=%u for %s" %
                (name_a, len(ra), name_b, len(rb), lab))
        for i, ((aa, xa, ya, za), (ba, xb, yb, zb)) in enumerate(zip(ra, rb)):
            if aa != ba:
                err("%s residue %u aa %s vs %s" % (lab, i, aa, ba))
            if (abs(xa - xb) > CCoordTol or abs(ya - yb) > CCoordTol
                    or abs(za - zb) > CCoordTol):
                err("%s residue %u coords %s vs %s" %
                    (lab, i, (xa, ya, za), (xb, yb, zb)))


def run_reseek(reseek, args):
    cmd = [reseek] + args
    print("+", " ".join(cmd))
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True)
    sys.stdout.write(r.stdout)
    if r.returncode != 0:
        err("reseek failed: %s" % cmd)


def main():
    reseek = os.environ.get("reseek", "../bin/reseek")
    if os.name == "nt" and not os.path.isfile(reseek):
        for cand in (
                "../src/Release/reseek.exe",
                "../Release/reseek.exe",
                "reseek.exe",
                ):
            if os.path.isfile(cand):
                reseek = cand
                break
    if not os.path.isfile(reseek) and not os.path.isfile(reseek + ".exe"):
        err("reseek not found: %s (set reseek=PATH)" % reseek)

    chains = [
        make_chain("chainA", 40, 0.0),
        make_chain("chainB", 200, 100.0),
        make_chain("chainC", 50, 200.0),
        make_chain("chainD", 80, 300.0),
        ]
    src_dict = dict(chains)

    tmp = tempfile.mkdtemp(prefix="cal_to_bcb_")
    src = os.path.join(tmp, "in.cal")
    write_cal(src, chains)

    bcb1 = os.path.join(tmp, "t1.bcb")
    bcb4 = os.path.join(tmp, "t4.bcb")
    cal1 = os.path.join(tmp, "from1.cal")
    cal4 = os.path.join(tmp, "from4.cal")

    run_reseek(reseek, [
        "-convert", src, "-bcb", bcb1,
        "-threads", "1", "-minchainlength", "1",
        "-log", os.path.join(tmp, "t1.log"),
        ])
    run_reseek(reseek, [
        "-convert", src, "-bcb", bcb4,
        "-threads", "4", "-minchainlength", "1",
        "-log", os.path.join(tmp, "t4.log"),
        ])
    run_reseek(reseek, [
        "-convert", bcb1, "-cal", cal1,
        "-threads", "1", "-minchainlength", "1",
        "-log", os.path.join(tmp, "from1.log"),
        ])
    run_reseek(reseek, [
        "-convert", bcb4, "-cal", cal4,
        "-threads", "1", "-minchainlength", "1",
        "-log", os.path.join(tmp, "from4.log"),
        ])

    d1 = read_cal_chains(cal1)
    d4 = read_cal_chains(cal4)
    compare_cal(src_dict, d1, "input", "1-thread roundtrip")
    compare_cal(src_dict, d4, "input", "4-thread roundtrip")
    compare_cal(d1, d4, "1-thread", "4-thread")
    print("OK %u chains  %s" % (len(src_dict), tmp))


if __name__ == "__main__":
    main()
