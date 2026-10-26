#!/usr/bin/env python3
"""Make a degraded annotation for the recovery benchmark (DESIGN.md §11.3).

Truth: the protein-coding genes of a GENCODE GTF, all their transcripts.
Degraded annotation: one transcript per gene (the Ensembl canonical coding one), with
known damage applied to a random subset of genes:

  deleted     the gene is removed                      (tests novel-gene recovery)
  utr_trim    the model is cut back to its CDS          (tests UTR recovery)
  split       the model is cut into two genes at a      (tests merging)
              central intron, each keeping its part of the CDS
  (all)       alternative isoforms are dropped          (tests new-isoform recovery)

Writes <out>.degraded.gtf, <out>.truth.gtf and <out>.damage.tsv (gene, damage).

usage: bench_degrade.py annotation out_prefix [--canonical gencode|first|longest_cds]
                        [--seed 1] [--deleted 0.1] [--utr_trim 0.2] [--split 0.05]

--canonical gencode (default): GENCODE GTF, protein_coding genes, Ensembl_canonical model.
--canonical first: any GTF/GFF3 (e.g. Araport11), genes with a coding model; canonical =
    the first coding model by ID (AT1G01010.1-style representative models).
--canonical longest_cds: as `first`, canonical = the model with the longest CDS.
"""

import argparse
import random
import re
from collections import defaultdict

from pasa_xg.gff import read_annotation, write_gtf
from pasa_xg.models import AnnotationSet, Gene, Transcript


def parse(gtf):
    tx = {}
    for line in open(gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or 'gene_type "protein_coding"' not in f[8]:
            continue
        tid = re.search(r'transcript_id "([^"]+)"', f[8])
        if not tid:
            continue
        tid = tid.group(1)
        r = tx.setdefault(tid, {"gene": re.search(r'gene_id "([^"]+)"', f[8]).group(1),
                                "contig": f[0], "strand": f[6], "exons": [], "cds": [],
                                "phase": {}, "canonical": False, "name": ""})
        if f[2] == "transcript":
            r["canonical"] = 'tag "Ensembl_canonical"' in f[8]
            m = re.search(r'gene_name "([^"]+)"', f[8])
            r["name"] = m.group(1) if m else ""
        elif f[2] == "exon":
            r["exons"].append((int(f[3]), int(f[4])))
        elif f[2] in ("CDS", "stop_codon"):
            r["cds"].append((int(f[3]), int(f[4])))
            if f[2] == "CDS":
                r["phase"][(int(f[3]), int(f[4]))] = int(f[7])
    out = defaultdict(list)
    for tid, r in tx.items():
        cds = _join(r["cds"])
        off = 0
        if r["phase"]:
            first = (min if r["strand"] == "+" else max)(r["phase"], key=lambda s: s[0])
            off = r["phase"][first]
        t = Transcript(tid, r["gene"], r["contig"], r["strand"], r["exons"], cds,
                       {"canonical": str(r["canonical"])}, cds_offset=off)
        out[r["gene"]].append(t)
    return out


def parse_generic(path, how):
    """Genes with at least one coding model; all their coding models, one marked canonical."""
    out = defaultdict(list)
    for g in read_annotation(path):
        coding = sorted((t for t in g.transcripts.values() if t.cds), key=lambda t: t.id)
        if not coding:
            continue
        canon = coding[0] if how == "first" else max(coding, key=lambda t: (t.cds_length(), t.id))
        for t in coding:
            out[g.id].append(t.copy(gene_id=g.id, attrs={"canonical": str(t is canon)}))
    return out


def _join(spans):
    o = []
    for a, b in sorted(spans):
        if o and a <= o[-1][1] + 1:
            o[-1][1] = max(o[-1][1], b)
        else:
            o.append([a, b])
    return [(a, b) for a, b in o]


def cds_only(t):
    lo, hi = t.cds[0][0], t.cds[-1][1]
    ex = [(max(a, lo), min(b, hi)) for a, b in t.exons if b >= lo and a <= hi]
    return t.copy(exons=ex)


def split(t, gid):
    """Two models from t at its central intron (each with its part of the CDS)."""
    k = len(t.exons) // 2
    left, right = t.exons[:k], t.exons[k:]
    lc = [s for s in t.cds if s[1] <= left[-1][1]]
    rc = [s for s in t.cds if s[0] >= right[0][0]]
    if not lc or not rc:
        return None
    a = Transcript(t.id + "_a", gid + "_a", t.contig, t.strand, left, lc)
    b = Transcript(t.id + "_b", gid + "_b", t.contig, t.strand, right, rc)
    # the downstream fragment (in transcription order) starts mid-codon
    up, down = (a, b) if t.strand == "+" else (b, a)
    up.cds_offset = t.cds_offset
    done = sum(e - s + 1 for s, e in up.cds) - up.cds_offset
    down.cds_offset = (3 - done % 3) % 3
    return a, b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gtf")
    ap.add_argument("out")
    ap.add_argument("--canonical", choices=("gencode", "first", "longest_cds"), default="gencode")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--deleted", type=float, default=0.10)
    ap.add_argument("--utr_trim", type=float, default=0.20)
    ap.add_argument("--split", type=float, default=0.05)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    genes = parse(a.gtf) if a.canonical == "gencode" else parse_generic(a.gtf, a.canonical)
    truth, degraded, damage = [], [], []
    for gid in sorted(genes):
        txs = genes[gid]
        g = Gene(gid, txs[0].contig, txs[0].strand)
        for t in txs:
            g.add(t.copy())
        truth.append(g)
        canon = [t for t in txs if t.attrs["canonical"] == "True" and t.cds]
        if not canon:
            continue
        c = canon[0].copy(attrs={})
        r = rng.random()
        if r < a.deleted:
            damage.append((gid, "deleted", c.id))
            continue
        r -= a.deleted
        if r < a.utr_trim and c.utrs() != ([], []):
            m = cds_only(c)
            damage.append((gid, "utr_trim", c.id))
            degraded.append(Gene(gid, c.contig, c.strand, {m.id: m}))
            continue
        r -= a.utr_trim
        if r < a.split and c.num_exons >= 4:
            parts = split(c, gid)
            if parts:
                damage.append((gid, "split", c.id))
                for p in parts:
                    degraded.append(Gene(p.gene_id, p.contig, p.strand, {p.id: p}))
                continue
        damage.append((gid, "none" if len(txs) == 1 else "isoforms_dropped", c.id))
        degraded.append(Gene(gid, c.contig, c.strand, {c.id: c}))
    with open(a.out + ".truth.gtf", "w") as fh:
        write_gtf(AnnotationSet(truth), fh)
    with open(a.out + ".degraded.gtf", "w") as fh:
        write_gtf(AnnotationSet(degraded), fh)
    with open(a.out + ".damage.tsv", "w") as fh:
        fh.write("gene_id\tdamage\tcanonical\n")
        for row in damage:
            fh.write("\t".join(row) + "\n")
    from collections import Counter
    print(Counter(d for _, d, _ in damage))


if __name__ == "__main__":
    main()
