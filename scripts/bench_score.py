#!/usr/bin/env python3
"""Score an updated annotation against the recovery benchmark (bench_degrade.py).

For each kind of damage, how much the update restored, next to how much the evidence
(the LRAA isoforms) could have restored ("supported"):

  deleted           novel gene whose intron chain equals a truth transcript of the gene
  utr_trim          model regains UTR on both ends with the truth intron chain
  split             the two fragments are joined into one gene with the truth chain
  isoforms_dropped  truth alternative isoforms (intron chains) restored

plus the precision of added or changed transcripts (intron chain equals any transcript
in the full GENCODE annotation of the region) and how many untouched canonical models
were changed away from the truth.

usage: bench_score.py prefix full_gencode.gtf updated.gff3 lraa_isoforms.gtf
"""

import sys
from collections import Counter, defaultdict

from pasa_xg.gff import read_annotation


def chain(t):
    return (t.contig, t.strand, tuple(t.introns())) if t.num_exons > 1 else None


def main(prefix, gencode, updated_p, lraa_p):
    truth = read_annotation(prefix + ".truth.gtf")
    full = read_annotation(gencode)
    degraded = read_annotation(prefix + ".degraded.gtf")
    upd = read_annotation(updated_p)
    lraa = read_annotation(lraa_p)
    damage = {}
    for line in open(prefix + ".damage.tsv").read().splitlines()[1:]:
        gid, kind, canon = line.split("\t")
        damage[gid] = (kind, canon)
    supported = {chain(t) for t in lraa.transcripts()} - {None}
    upd_chains = defaultdict(set)
    for g in upd:
        for t in g.transcripts.values():
            upd_chains[chain(t)].add(g.id)
    all_truth = {chain(t) for t in full.transcripts()} - {None}
    res = Counter()

    for gid, (kind, canon) in damage.items():
        tg = truth.genes[gid]
        chains = {chain(t) for t in tg.transcripts.values()} - {None}
        ct = tg.transcripts[canon]
        if kind == "deleted":
            res["deleted"] += 1
            if chains & supported:
                res["deleted.supported"] += 1
            if any(c in upd_chains and any(h not in degraded.genes for h in upd_chains[c]) for c in chains):
                res["deleted.recovered"] += 1
        elif kind == "utr_trim":
            res["utr_trim"] += 1
            if chain(ct) in supported or ct.num_exons == 1:
                res["utr_trim.supported"] += 1
            u = upd.genes.get(gid)
            m = u.transcripts.get(canon) if u else None
            if m is not None and m.cds:
                five, three = m.utrs()
                if five and three:
                    res["utr_trim.both_utrs"] += 1
                if (five or three):
                    res["utr_trim.any_utr"] += 1
        elif kind == "split":
            res["split"] += 1
            if chain(ct) in supported:
                res["split.supported"] += 1
            if not (gid + "_a" in upd.genes and gid + "_b" in upd.genes):
                res["split.joined"] += 1
                if chain(ct) in upd_chains:
                    res["split.joined_truth_chain"] += 1
        if kind in ("isoforms_dropped", "utr_trim", "split", "none"):
            alts = {chain(t) for t in tg.transcripts.values() if t.id != canon} - {None, chain(ct)}
            res["alt_isoforms"] += len(alts)
            res["alt_isoforms.supported"] += len(alts & supported)
            res["alt_isoforms.recovered"] += len(alts & set(upd_chains))
        if kind in ("isoforms_dropped", "none"):
            res["untouched_models"] += 1
            u = upd.genes.get(gid)
            m = u.transcripts.get(canon) if u else None
            if m is None:
                res["untouched_models.removed"] += 1
            elif chain(m) != chain(ct):
                res["untouched_models.chain_changed"] += 1
            elif m.cds != ct.cds:
                res["untouched_models.cds_changed"] += 1

    deg_keys = {(t.id, chain(t), tuple(t.cds)) for t in degraded.transcripts()}
    added = [t for t in upd.transcripts() if (t.id, chain(t), tuple(t.cds)) not in deg_keys]
    multi = [t for t in added if chain(t)]
    res["added_or_changed"] = len(added)
    res["added_or_changed.multi_exon"] = len(multi)
    res["added_or_changed.multi_exon.in_gencode"] = sum(1 for t in multi if chain(t) in all_truth)

    def pct(a, b):
        return f"{100.0 * res[a] / res[b]:.0f}%" if res[b] else "-"

    print(f"{'measure':42s} {'count':>7s}  of")
    for k in sorted(res):
        base = k.rsplit(".", 1)[0] if "." in k else None
        print(f"{k:42s} {res[k]:7d}  {pct(k, base) if base else ''}")


if __name__ == "__main__":
    main(*sys.argv[1:])
