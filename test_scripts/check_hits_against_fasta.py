#!/usr/bin/env python3

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fasta import ReadSeqsDict

PCTID_TOL = 0.05

VALID_COLUMNS = frozenset([
    "query", "target", "pvalue", "evalue",
    "qlo", "qhi", "tlo", "thi", "ql", "tl",
    "pctid", "cigar", "qrow", "trow", "qrowg", "trowg",
    "raw", "mega", "lddt", "dali", "ids", "gaps", "aq", "cols",
    "qcovpct", "tcovpct",
])

errors = 0


def err(msg):
    global errors
    errors += 1
    print("ERROR: %s" % msg)


def is_upper_aa(c):
    return "A" <= c <= "Z"


def is_lower_aa(c):
    return "a" <= c <= "z"


def column_kind(q, t):
    if is_upper_aa(q) and is_upper_aa(t):
        return "aligned"
    if q == "-" and is_upper_aa(t):
        return "aligned"
    if t == "-" and is_upper_aa(q):
        return "aligned"
    if (is_lower_aa(q) or q == ".") and (is_lower_aa(t) or t == "."):
        return "unaligned"
    return "invalid"


def strip_gaps(row):
    return "".join(c for c in row if c not in ".-")


def parse_cigar(cigar):
    if not cigar:
        return 0, "", False
    ns = []
    ops = []
    n = 0
    for c in cigar:
        if c.isdigit():
            n = n * 10 + int(c)
        else:
            ops.append(c)
            ns.append(n)
            n = 0
    if n != 0:
        raise ValueError("trailing digits in cigar '%s'" % cigar)
    if len(ops) != len(ns):
        raise ValueError("bad cigar '%s'" % cigar)
    soft_clip_q = 0
    has_s = False
    path = []
    for op, cnt in zip(ops, ns):
        if op == "S":
            soft_clip_q = cnt
            has_s = True
        elif op == "T":
            pass
        elif op in "MDI":
            path.extend([op] * cnt)
        else:
            raise ValueError("unsupported cigar op '%s' in '%s'" % (op, cigar))
    return soft_clip_q, "".join(path), has_s


def rows_from_path(path, qseq, tseq, qlo0, tlo0):
    qpos = qlo0
    tpos = tlo0
    qrow = []
    trow = []
    for c in path:
        if c == "M":
            if qpos >= len(qseq) or tpos >= len(tseq):
                raise ValueError("cigar path exceeds sequence bounds")
            qrow.append(qseq[qpos])
            trow.append(tseq[tpos])
            qpos += 1
            tpos += 1
        elif c == "D":
            if tpos >= len(tseq):
                raise ValueError("cigar D exceeds target bounds")
            qrow.append("-")
            trow.append(tseq[tpos])
            tpos += 1
        elif c == "I":
            if qpos >= len(qseq):
                raise ValueError("cigar I exceeds query bounds")
            qrow.append(qseq[qpos])
            trow.append("-")
            qpos += 1
        else:
            raise ValueError("bad path op '%s'" % c)
    return "".join(qrow), "".join(trow), qpos - 1, tpos - 1


def count_ids_gaps(qrow, trow):
    ids = 0
    gaps = 0
    aligned_cols = 0
    for q, t in zip(qrow, trow):
        kind = column_kind(q, t)
        if kind != "aligned":
            continue
        aligned_cols += 1
        if q == t and is_upper_aa(q):
            ids += 1
        if (q == "-") ^ (t == "-"):
            gaps += 1
    return ids, gaps, aligned_cols


def validate_row_chars(row, allow_unaligned, row_name, lineno):
    for c in row:
        if is_upper_aa(c) or c == "-":
            continue
        if allow_unaligned and (is_lower_aa(c) or c == "."):
            continue
        err("line %u: %s has invalid char '%s'" % (lineno, row_name, c))
        return False
    return True


def validate_row_pair(qrow, trow, allow_unaligned, lineno):
    if len(qrow) != len(trow):
        err("line %u: qrow len %u != trow len %u" % (lineno, len(qrow), len(trow)))
        return None
    if not validate_row_chars(qrow, allow_unaligned, "qrow", lineno):
        return None
    if not validate_row_chars(trow, allow_unaligned, "trow", lineno):
        return None
    for i, (q, t) in enumerate(zip(qrow, trow)):
        kind = column_kind(q, t)
        if kind == "invalid":
            err("line %u: col %u invalid pair (%s, %s)" % (lineno, i, q, t))
            return None
        if not allow_unaligned and kind == "unaligned":
            err("line %u: col %u unaligned not allowed in qrow/trow" % (lineno, i))
            return None
    return count_ids_gaps(qrow, trow)


def validate_row_substring(row, seq, lo, hi, row_name, lineno):
    if lo is None or hi is None:
        err("line %u: %s present but qlo/qhi or tlo/thi missing" % (lineno, row_name))
        return
    seg = strip_gaps(row)
    expected = seq[lo - 1:hi]
    if seg != expected:
        err("line %u: %s segment mismatch (got len %u, expected len %u)" % (
            lineno, row_name, len(seg), len(expected)))


def validate_coords(lo, hi, length, lo_name, hi_name, length_name, lineno):
    if lo is None or hi is None or length is None:
        return
    if lo < 1 or hi < lo or hi > length:
        err("line %u: bad %s=%s %s=%s %s=%u" % (
            lineno, lo_name, lo, hi_name, hi, length_name, length))


def validate_cigar_rows(cigar, qrow, trow, qseq, tseq, qlo, tlo, lineno):
    try:
        soft_clip_q, path, has_s = parse_cigar(cigar)
    except ValueError as e:
        err("line %u: %s" % (lineno, e))
        return
    if has_s and qlo is not None and soft_clip_q != qlo - 1:
        err("line %u: cigar S=%u but qlo=%u (expected S=%u)" % (
            lineno, soft_clip_q, qlo, qlo - 1))
    if len(path) != len(qrow):
        err("line %u: cigar path len %u != row len %u" % (lineno, len(path), len(qrow)))
        return
    qlo0 = qlo - 1 if qlo is not None else soft_clip_q
    tlo0 = tlo - 1 if tlo is not None else 0
    try:
        exp_qrow, exp_trow, qhi0, thi0 = rows_from_path(path, qseq, tseq, qlo0, tlo0)
    except ValueError as e:
        err("line %u: %s" % (lineno, e))
        return
    if exp_qrow != qrow:
        err("line %u: cigar/qrow mismatch" % lineno)
    if exp_trow != trow:
        err("line %u: cigar/trow mismatch" % lineno)
    return qhi0 + 1, thi0 + 1


def parse_int(val, name, lineno):
    try:
        return int(val)
    except ValueError:
        err("line %u: %s '%s' is not an integer" % (lineno, name, val))
        return None


def parse_float(val, name, lineno):
    try:
        return float(val)
    except ValueError:
        err("line %u: %s '%s' is not a float" % (lineno, name, val))
        return None


def validate_line(fields, col_idx, seqs, lineno):
    def col(name):
        if name not in col_idx:
            return None
        return fields[col_idx[name]]

    query = col("query")
    target = col("target")
    qseq = None
    tseq = None

    if query is not None:
        if query not in seqs:
            err("line %u: query '%s' not in FASTA" % (lineno, query))
        else:
            qseq = seqs[query]
    if target is not None:
        if target not in seqs:
            err("line %u: target '%s' not in FASTA" % (lineno, target))
        else:
            tseq = seqs[target]

    ql = parse_int(col("ql"), "ql", lineno) if "ql" in col_idx else None
    tl = parse_int(col("tl"), "tl", lineno) if "tl" in col_idx else None
    qlo = parse_int(col("qlo"), "qlo", lineno) if "qlo" in col_idx else None
    qhi = parse_int(col("qhi"), "qhi", lineno) if "qhi" in col_idx else None
    tlo = parse_int(col("tlo"), "tlo", lineno) if "tlo" in col_idx else None
    thi = parse_int(col("thi"), "thi", lineno) if "thi" in col_idx else None

    if ql is not None and qseq is not None and ql != len(qseq):
        err("line %u: ql=%u != len(query)=%u" % (lineno, ql, len(qseq)))
    if tl is not None and tseq is not None and tl != len(tseq):
        err("line %u: tl=%u != len(target)=%u" % (lineno, tl, len(tseq)))

    validate_coords(qlo, qhi, ql, "qlo", "qhi", "ql", lineno)
    validate_coords(tlo, thi, tl, "tlo", "thi", "tl", lineno)

    qrow = col("qrow")
    trow = col("trow")
    qrowg = col("qrowg")
    trowg = col("trowg")

    if qrow is not None and qseq is not None:
        validate_row_substring(qrow, qseq, qlo, qhi, "qrow", lineno)
    if trow is not None and tseq is not None:
        validate_row_substring(trow, tseq, tlo, thi, "trow", lineno)

    if qrowg is not None and qseq is not None:
        if strip_gaps(qrowg).upper() != qseq.upper() and strip_gaps(qrowg) != qseq:
            seg = strip_gaps(qrowg)
            if seg != qseq:
                err("line %u: qrowg does not match query sequence" % lineno)
    if trowg is not None and tseq is not None:
        seg = strip_gaps(trowg)
        if seg != tseq:
            err("line %u: trowg does not match target sequence" % lineno)

    stats = None
    if qrow is not None and trow is not None:
        stats = validate_row_pair(qrow, trow, allow_unaligned=False, lineno=lineno)
    elif qrowg is not None and trowg is not None:
        stats = validate_row_pair(qrowg, trowg, allow_unaligned=True, lineno=lineno)

    if stats is not None:
        exp_ids, exp_gaps, aligned_cols = stats
        if "ids" in col_idx:
            ids = parse_int(col("ids"), "ids", lineno)
            if ids is not None and ids != exp_ids:
                err("line %u: ids=%u expected %u" % (lineno, ids, exp_ids))
        if "gaps" in col_idx:
            gaps = parse_int(col("gaps"), "gaps", lineno)
            if gaps is not None and gaps != exp_gaps:
                err("line %u: gaps=%u expected %u" % (lineno, gaps, exp_gaps))
        if "pctid" in col_idx:
            pctid = parse_float(col("pctid"), "pctid", lineno)
            if pctid is not None and aligned_cols > 0:
                exp_pctid = 100.0 * exp_ids / aligned_cols
                if abs(pctid - exp_pctid) > PCTID_TOL:
                    err("line %u: pctid=%.4f expected %.4f" % (lineno, pctid, exp_pctid))
        if "cols" in col_idx:
            cols = parse_int(col("cols"), "cols", lineno)
            row = qrow if qrow is not None else qrowg
            if cols is not None and row is not None and cols != len(row):
                err("line %u: cols=%u expected %u" % (lineno, cols, len(row)))

    cigar = col("cigar")
    if cigar is not None and qrow is not None and trow is not None:
        if qseq is None or tseq is None:
            err("line %u: cigar+qrow+trow need query and target in FASTA" % lineno)
        else:
            inferred = validate_cigar_rows(
                cigar, qrow, trow, qseq, tseq, qlo, tlo, lineno)
            if inferred is not None:
                qhi_inf, thi_inf = inferred
                if qhi is not None and qhi != qhi_inf:
                    err("line %u: qhi=%u expected %u from cigar" % (lineno, qhi, qhi_inf))
                if thi is not None and thi != thi_inf:
                    err("line %u: thi=%u expected %u from cigar" % (lineno, thi, thi_inf))


def usage():
    print("Usage: %s <fasta> <columns> <tsv>" % sys.argv[0], file=sys.stderr)
    print("  columns: + -separated names from userfieldnames.h", file=sys.stderr)
    print("  e.g. query+target+qlo+qhi+ql+tlo+thi+tl+cigar+qrow+trow+ids+gaps+pctid",
          file=sys.stderr)


def main():
    if len(sys.argv) != 4:
        usage()
        sys.exit(1)

    fasta_fn = sys.argv[1]
    columns = sys.argv[2].split("+")
    tsv_fn = sys.argv[3]

    if not columns or columns == [""]:
        err("empty columns list")
        sys.exit(1)

    for name in columns:
        if name not in VALID_COLUMNS:
            err("unknown column '%s'" % name)
    if errors:
        sys.exit(1)

    col_idx = {name: i for i, name in enumerate(columns)}
    expected_nflds = len(columns)

    if not os.path.isfile(fasta_fn):
        err("FASTA not found: %s" % fasta_fn)
        sys.exit(1)
    if not os.path.isfile(tsv_fn):
        err("TSV not found: %s" % tsv_fn)
        sys.exit(1)

    seqs = ReadSeqsDict(fasta_fn)

    nlines = 0
    for lineno, line in enumerate(open(tsv_fn), 1):
        line = line.rstrip("\n\r")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != expected_nflds:
            err("line %u: expected %u fields, got %u" % (lineno, expected_nflds, len(fields)))
            continue
        nlines += 1
        validate_line(fields, col_idx, seqs, lineno)

    if nlines == 0:
        err("no data lines in %s" % tsv_fn)

    if errors:
        print("%s FAILED: %u error(s)" % (sys.argv[0], errors), file=sys.stderr)
        sys.exit(1)

    print("%s SUCCESS (%u lines)" % (sys.argv[0], nlines))


if __name__ == "__main__":
    main()
