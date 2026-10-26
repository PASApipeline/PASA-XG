#!/usr/bin/env python3
"""Side-by-side table of bench_score.py outputs.

usage: bench_table.py name=score.txt [name=score.txt ...]
"""

import sys

ROWS = [("deleted genes recovered", "deleted.recovered", "deleted.supported"),
        ("UTRs restored (any)", "utr_trim.any_utr", "utr_trim.supported"),
        ("split genes rejoined, true chain", "split.joined_truth_chain", "split.supported"),
        ("alt isoforms restored", "alt_isoforms.recovered", "alt_isoforms.supported"),
        ("correct models: chain changed", "untouched_models.chain_changed", "untouched_models"),
        ("correct models: CDS changed", "untouched_models.cds_changed", "untouched_models"),
        ("added/changed multi-exon", "added_or_changed.multi_exon", None),
        ("  of which in GENCODE", "added_or_changed.multi_exon.in_gencode", None)]


def load(p):
    out = {}
    for line in open(p).read().splitlines()[1:]:
        f = line.split()
        out[f[0]] = int(f[1])
    return out


cols = [(a.split("=", 1)[0], load(a.split("=", 1)[1])) for a in sys.argv[1:]]
w = max(len(n) for n, _ in cols) + 2
print(f"{'':36s}" + "".join(f"{n:>{w}s}" for n, _ in cols) + "   of")
for label, k, base in ROWS:
    ref = cols[0][1].get(base, "") if base else ""
    print(f"{label:36s}" + "".join(f"{c.get(k, 0):>{w}d}" for _, c in cols) + f"   {ref}")
