#!/usr/bin/env python3
"""Compare a PASA-XG updated annotation with a Perl PASA update of the same input.

Reports, against the original annotation:
  - which genes each pipeline changed (structure, isoforms added, merged/split),
  - transcript structures (intron chain + CDS) found by both, by one only,
  - novel genes from each, and how many overlap each other.

usage: compare_with_perl.py original.gff3 perl_updated.gff3 pasa_xg_updated.gff3
"""

import sys
from collections import Counter

from pasa_xg.gff import read_annotation


def key(t, with_cds=True):
    return (t.contig, t.strand, tuple(t.introns()) if t.num_exons > 1 else (t.start, t.end),
            tuple(t.cds) if with_cds else ())


def structures(ann, with_cds=True):
    return {key(t, with_cds): t.id for t in ann.transcripts()}


def main(orig_p, perl_p, xg_p):
    orig, perl, xg = (read_annotation(p) for p in (orig_p, perl_p, xg_p))
    print(f"genes: original {len(orig)}, Perl {len(perl)}, PASA-XG {len(xg)}")
    print(f"transcripts: original {sum(1 for _ in orig.transcripts())}, "
          f"Perl {sum(1 for _ in perl.transcripts())}, PASA-XG {sum(1 for _ in xg.transcripts())}")
    o = structures(orig)
    for label, with_cds in (("intron chain + CDS", True), ("intron chain only", False)):
        p, x = structures(perl, with_cds), structures(xg, with_cds)
        o2 = structures(orig, with_cds)
        new_p, new_x = set(p) - set(o2), set(x) - set(o2)
        print(f"\nnew or changed structures ({label}): Perl {len(new_p)}, PASA-XG {len(new_x)}, "
              f"both {len(new_p & new_x)}")
        kept_p = len(set(o2) & set(p))
        kept_x = len(set(o2) & set(x))
        print(f"original structures kept unchanged: Perl {kept_p}, PASA-XG {kept_x} (of {len(o2)})")

    def changed_genes(upd):
        out = Counter()
        for g in orig:
            if g.id not in upd.genes:
                out["removed (merged)"] += 1
                continue
            u = upd.genes[g.id]
            ok = {key(t) for t in g.transcripts.values()}
            uk = {key(t) for t in u.transcripts.values()}
            if ok == uk:
                out["unchanged"] += 1
            elif len(u.transcripts) > len(g.transcripts):
                out["isoforms added"] += 1
            else:
                out["structure updated"] += 1
        return out

    cp, cx = changed_genes(perl), changed_genes(xg)
    print("\noriginal genes:")
    for k in sorted(set(cp) | set(cx)):
        print(f"  {k:20s} Perl {cp[k]:5d}  PASA-XG {cx[k]:5d}")
    # gene-level agreement on which genes changed
    def changed_set(upd):
        s = set()
        for g in orig:
            u = upd.genes.get(g.id)
            if u is None or {key(t) for t in g.transcripts.values()} != {key(t) for t in u.transcripts.values()}:
                s.add(g.id)
        return s
    sp, sx = changed_set(perl), changed_set(xg)
    print(f"changed by both: {len(sp & sx)}, Perl only: {len(sp - sx)}, PASA-XG only: {len(sx - sp)}")
    nov_p = [g for g in perl if g.id not in orig.genes]
    nov_x = [g for g in xg if g.id not in orig.genes]
    ov = sum(1 for g in nov_x if any(h.strand == g.strand and h.start <= g.end and g.start <= h.end
                                     and h.contig == g.contig for h in nov_p))
    print(f"\nnew genes: Perl {len(nov_p)}, PASA-XG {len(nov_x)}; PASA-XG new genes overlapping a Perl new gene: {ov}")


if __name__ == "__main__":
    main(*sys.argv[1:])
