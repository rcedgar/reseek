#!/usr/bin/env python3

import os
import subprocess
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.normpath(os.path.join(script_dir, ".."))
test_data = os.path.join(repo_root, "test_data")
outdir = os.path.join(repo_root, "test_output", "label_by_filename")

reseek = os.environ.get("reseek")
if not reseek:
    for candidate in (
        os.path.join(repo_root, "src", "Release", "reseek.exe"),
        os.path.join(repo_root, "src", "Release", "reseek"),
        os.path.join(repo_root, "bin", "reseek.exe"),
        os.path.join(repo_root, "bin", "reseek"),
    ):
        if os.path.isfile(candidate):
            reseek = candidate
            break

errors = 0


def err(msg):
    global errors
    errors += 1
    print("ERROR: %s" % msg)


def read_fasta_labels(fn):
    labels = []
    for line in open(fn):
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            labels.append(line[1:])
    return labels


def run_convert(struct_path, extra_args, out_fa):
    cmd = [reseek, "-convert", struct_path, "-minchainlength", "3"] + extra_args + ["-fasta", out_fa]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        err("command failed (%d): %s\n%s" % (
            result.returncode, " ".join(cmd), result.stderr))
    return result


def check_label_prefix(labels, expected_prefix, context):
    if not labels:
        err("%s: no FASTA labels" % context)
        return
    label = labels[0]
    if not label.startswith(expected_prefix):
        err("%s: expected label prefix '%s', got '%s'" % (
            context, expected_prefix, label))


def main():
    global reseek
    if not reseek or not os.path.isfile(reseek):
        err("reseek binary not found; set reseek=PATH or build first")
        print("%s FAILED: %u error(s)" % (sys.argv[0], errors))
        sys.exit(1)

    os.makedirs(outdir, exist_ok=True)

    pdb_path = os.path.join(test_data, "label_by_fn.pdb")
    cif_path = os.path.join(test_data, "label_by_fn.cif")

    for fn in (pdb_path, cif_path):
        if not os.path.isfile(fn):
            err("missing test file: %s" % fn)

    if errors:
        print("%s FAILED: %u error(s)" % (sys.argv[0], errors))
        sys.exit(1)

    pdb_default = os.path.join(outdir, "pdb_default.fa")
    pdb_by_fn = os.path.join(outdir, "pdb_by_fn.fa")
    cif_default = os.path.join(outdir, "cif_default.fa")
    cif_by_fn = os.path.join(outdir, "cif_by_fn.fa")

    run_convert(pdb_path, [], pdb_default)
    run_convert(pdb_path, ["-label_by_filename"], pdb_by_fn)
    run_convert(cif_path, [], cif_default)
    run_convert(cif_path, ["-label_by_filename"], cif_by_fn)

    check_label_prefix(read_fasta_labels(pdb_default), "ABCD_A", "PDB default")
    check_label_prefix(read_fasta_labels(pdb_by_fn), "label_by_fn_A", "PDB -label_by_filename")
    check_label_prefix(read_fasta_labels(cif_default), "TESTID_A", "CIF default")
    check_label_prefix(read_fasta_labels(cif_by_fn), "label_by_fn_A", "CIF -label_by_filename")

    if errors:
        print("%s FAILED: %u error(s)" % (sys.argv[0], errors))
        sys.exit(1)

    print("%s SUCCESS" % sys.argv[0])


if __name__ == "__main__":
    main()
