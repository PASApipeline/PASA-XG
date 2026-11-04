"""Noncoding pass (--include_noncoding, docs/NONCODING_DESIGN.md).

Runs after the coding rounds on the final annotation. It takes the isoforms the coding rounds did
not incorporate and adds those with a new structure and enough evidence:

- an isoform that belongs to exactly one gene on its strand (assign_genes: shares an intron
  with it, or else overlaps its exons, small-RNA genes aside) becomes a transcript of that gene, unless it assembles with one of the gene's transcripts (merge_compatible): then it is a
  match, fragment or end variant of that transcript, and nothing is added;
- a spliced isoform overlapping no gene on its strand joins a novel noncoding gene (isoforms that
  overlap each other share one);
- an isoform overlapping two or more genes on its strand (readthrough) is left out.

Labels: pasa_xg_coding "low" (the isoform's ORF is the gene's product or a contiguous part of it,
and not an NMD target: the CDS is kept), "nmd" (complete ORF, stop codon upstream of the last
junction) or "none". Every gene gets pasa_xg_gene_class "coding" or "likely_noncoding".

Evidence: n_samples from the isoform's GTF attribute (missing: 1, as for a single-sample run) and
supported ends (an isoform TSS / PolyA flag, or a known site within the window).
  low                    n_samples >= 2, or >= 1 supported end
  other spliced          n_samples >= 2 and >= 1 supported end, or n_samples >= 3
  mono-exonic            n_samples >= 3, both ends supported, not inside a reference intron
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

from ..genome import Genome
from ..isoforms import Isoform
from ..models import AnnotationSet, Gene, Transcript
from . import checks as C
from .apply import Namer
from .sites import POLYA, TSS, SiteIndex, five_prime, three_prime
from .structure import introns_of, merge_compatible

log = logging.getLogger("pasa_xg")

# decision statuses of the noncoding pass (beside PASA's, compare.STATUS)
NC_STATUS = {
    101: "noncoding pass: isoform added to a known gene",
    102: "noncoding pass: isoform of a novel noncoding gene",
    103: "noncoding pass: not enough evidence",
    104: "noncoding pass: assembles with an existing transcript",
    105: "noncoding pass: overlaps two or more genes",
    106: "noncoding pass: mono-exonic isoform outside a known gene",
}
NC_ACTION = {101: "noncoding_isoform", 102: "noncoding_gene"}


@dataclass
class NcDecision:
    isoform_id: str
    status: int
    comment: str
    genes: List[str]
    annotation_id: str = ""

    @property
    def action(self) -> str:
        return NC_ACTION.get(self.status, "rejected")


def n_samples(iso: Isoform) -> int:
    try:
        return int(float(iso.tx.attrs.get("n_samples", 1)))
    except ValueError:
        return 1


def supported_ends(iso: Isoform, sites: Optional[SiteIndex]) -> Tuple[bool, bool]:
    t = iso.tx
    tss, pa = iso.tss, iso.polya
    if sites is not None:
        tss = tss or sites.near(t.contig, t.strand, TSS, five_prime(t)) is not None
        pa = pa or sites.near(t.contig, t.strand, POLYA, three_prime(t)) is not None
    return tss, pa


def evidence_ok(iso: Isoform, coding: str, sites: Optional[SiteIndex], in_intron: bool) -> Tuple[bool, str]:
    n = n_samples(iso)
    tss, pa = supported_ends(iso, sites)
    ends = tss + pa
    if iso.tx.num_exons == 1:
        ok = n >= 3 and ends == 2 and not in_intron
    elif coding == "low":
        ok = n >= 2 or ends >= 1
    else:
        ok = (n >= 2 and ends >= 1) or n >= 3
    why = f"{n} sample(s), {ends} supported end(s)" + (", inside a reference intron" if in_intron else "")
    return ok, why


def coding_label(iso: Isoform, models: List[Transcript], genome: Genome, table, p: C.Params
                 ) -> Tuple[str, Optional[Transcript]]:
    """('low', isoform with its CDS) when its ORF is the gene's product or a contiguous part of
    it and is not an NMD target; else 'nmd' or 'none' (and no CDS)."""
    if not iso.tx.cds:          # the coding stage puts the isoform's ORF on it as its CDS
        return "none", None
    t = iso.tx.copy()
    pep = genome.protein(t, table)
    nmd = pep.endswith("*") and C.nmd_distance(t) > p.NMD_STOP_TO_JUNCTION
    if nmd:
        return "nmd", None
    core = pep.rstrip("*")
    for m in models:
        ref = genome.protein(m, table).rstrip("*")
        if core and core in ref:
            return "low", t
        ok, _ = C.compare_proteins(m, t, genome, table, p, True, not pep.startswith("M"), not pep.endswith("*"))
        if ok:
            return "low", t
    return "none", None


def _exon_overlap(x: Iterable[Tuple[int, int]], y: Iterable[Tuple[int, int]]) -> bool:
    return any(a <= d and c <= b for a, b in x for c, d in y)


SMALL_RNA_LEN = 300   # genes whose transcripts are all shorter (miRNA, snoRNA, ...) are not hosts


def assign_genes(tx: Transcript, ann: AnnotationSet) -> List[Gene]:
    """Same-strand genes the isoform belongs to: of those sharing an intron with it, the ones whose
    best transcript's intron set is most similar (ties: all of them); if none shares an intron,
    those whose exons it overlaps, leaving out small-RNA genes embedded in a host gene."""
    over = [g for g in ann.overlapping(tx.contig, tx.start, tx.end, strand=tx.strand)
            if any(_exon_overlap(r.exons, tx.exons) for r in g.transcripts.values())]
    ints = set(introns_of(tx.exons))
    if ints:
        # intron-set similarity (Jaccard) of each gene's best transcript; overlapping genes that
        # share introns (a readthrough and its parts) go to the clearly better match
        best: Dict[str, float] = {}
        for g in ann.overlapping(tx.contig, tx.start, tx.end, strand=tx.strand):
            for r in g.transcripts.values():
                ri = set(introns_of(r.exons))
                if ints & ri:
                    best[g.id] = max(best.get(g.id, 0.0), len(ints & ri) / len(ints | ri))
        if best:
            top = max(best.values())
            return [ann.genes[k] for k, v in best.items() if v >= top - 1e-9]
    big = [g for g in over if max(r.length() for r in g.transcripts.values()) >= SMALL_RNA_LEN]
    return big if big or not over else over


def _in_intron(t: Transcript, ann: AnnotationSet) -> bool:
    lo, hi = t.start, t.end
    for g in ann.overlapping(t.contig, lo, hi):
        for r in g.transcripts.values():
            if any(a <= lo and hi <= b for a, b in introns_of(r.exons)):
                return True
    return False


def run_noncoding(ann: AnnotationSet, isos: Dict[str, Isoform], candidates: Iterable[str],
                  genome: Genome, table, p: C.Params, sites: Optional[SiteIndex], run_label: str = ""
                  ) -> Tuple[AnnotationSet, List[NcDecision]]:
    new = ann.copy(label="noncoding")
    namer = Namer(new)
    out: List[NcDecision] = []
    order = sorted(candidates, key=lambda k: (-n_samples(isos[k]), -isos[k].tpm, k))
    novel: List[Tuple[Transcript, Isoform]] = []

    def tag(t: Transcript, iso: Isoform, coding: str) -> None:
        tss, pa = supported_ends(iso, sites)
        t.attrs = {"transcript_type": "pasa_low_coding" if coding == "low" else "pasa_noncoding",
                   "pasa_xg_coding": coding, "pasa_xg_update": "noncoding_isoform",
                   "pasa_xg_round": "noncoding", "pasa_xg_isoform": iso.id,
                   **({"pasa_xg_run": run_label} if run_label else {}),
                   "n_samples": str(n_samples(iso)), "TPM": f"{iso.tpm:.3f}"}

    for k in order:
        iso = isos[k]
        tx = iso.tx
        genes = assign_genes(tx, new)
        if len(genes) > 1:
            out.append(NcDecision(k, 105, "", [g.id for g in genes]))
            continue
        if not genes:
            if tx.num_exons == 1:
                out.append(NcDecision(k, 106, "", []))
                continue
            label, _ = coding_label(iso, [], genome, table, p)
            ok, why = evidence_ok(iso, label, sites, False)
            if not ok:
                out.append(NcDecision(k, 103, why, []))
                continue
            hit = next((t for t, _ in novel if t.strand == tx.strand and t.contig == tx.contig
                        and merge_compatible(t.exons, tx.exons) is not None), None)
            if hit is not None:
                out.append(NcDecision(k, 104, f"assembles with {hit.id}", []))
                continue
            t = Transcript(k, "", tx.contig, tx.strand, list(tx.exons))
            tag(t, iso, label)
            novel.append((t, iso))
            out.append(NcDecision(k, 102, why, []))
            continue
        g = genes[0]
        hit = next((r for r in g.transcripts.values() if merge_compatible(r.exons, tx.exons) is not None), None)
        if hit is not None:
            out.append(NcDecision(k, 104, f"assembles with {hit.id}", [g.id]))
            continue
        label, coded = coding_label(iso, [r for r in g.transcripts.values() if r.cds], genome, table, p)
        ok, why = evidence_ok(iso, label, sites, tx.num_exons == 1 and _in_intron(tx, new))
        if not ok:
            out.append(NcDecision(k, 103, why, [g.id]))
            continue
        tid = namer.isoform(g.id)
        t = (coded or tx).copy(id=tid, gene_id=g.id, attrs={})
        if coded is None:
            t.cds = []
        tag(t, iso, label)
        if t.start < g.start or t.end > g.end:
            new._index = None   # the gene's span grows: rebuild the overlap index
        g.transcripts[tid] = t
        out.append(NcDecision(k, 101, f"{label}; {why}", [g.id], tid))

    # novel noncoding genes: isoforms sharing exonic overlap on a strand form one gene
    novel.sort(key=lambda x: (x[0].contig, x[0].strand, x[0].start))
    groups: List[List[Tuple[Transcript, Isoform]]] = []
    for t, iso in novel:
        if groups and any(u.contig == t.contig and u.strand == t.strand and _exon_overlap(u.exons, t.exons)
                          for u, _ in groups[-1]):
            groups[-1].append((t, iso))
        else:
            groups.append([(t, iso)])
    where: Dict[str, Tuple[str, str]] = {}
    for grp in groups:
        gid = namer.gene()
        tx = {}
        for i, (t, iso) in enumerate(grp, 1):
            t2 = t.copy(id=f"{gid}.t{i}", gene_id=gid)
            tx[t2.id] = t2
            where[iso.id] = (gid, t2.id)
        flip = "-" if grp[0][0].strand == "+" else "+"
        anti = sorted({g.id for t, _ in grp for g in ann.overlapping(t.contig, t.start, t.end) if g.strand == flip})
        new.add_gene(Gene(gid, grp[0][0].contig, grp[0][0].strand, tx,
                          {"pasa_xg_novel": "noncoding", **({"antisense_to": ",".join(anti)} if anti else {})}))
    for d in out:
        if d.status == 102:
            d.genes, d.annotation_id = [where[d.isoform_id][0]], where[d.isoform_id][1]
    for g in new:
        g.attrs = {**g.attrs, "pasa_xg_gene_class": "coding" if any(t.cds for t in g.transcripts.values())
                   else "likely_noncoding"}
    n = defaultdict(int)
    for d in out:
        n[d.action] += 1
    log.info("noncoding pass: %d candidates, %s", len(out), dict(n))
    return new, out
