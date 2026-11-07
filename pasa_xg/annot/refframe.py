"""Reference-frame check of ORFs, against the gene's reference coding transcripts.

pytransdecoder (and the update's own ORF calls on merged structures) take the best ORF of a sequence
on its own. Against the reference coding transcripts that share introns with it, two cases call for
a different reading:

- The transcript holds the best-matching reference's start codon upstream of its ORF. A scanning
  ribosome meets that ATG first, and here it opens a different, overlapping ORF. When the
  ORF from the reference ATG is long enough and not an NMD target it becomes the ORF; otherwise the
  transcript is taken as noncoding (no CDS).
- The transcript is end-incomplete in the reference CDS (best-matching reference): it starts after the
  reference start codon with no TSS support, or ends before the reference stop codon, with no polyA
  support and no stop codon in the reference frame. Its ORF is then read in the reference frame (from
  the first in-frame codon) to the first stop, so it is partial (cds_start_NF / cds_end_NF), and a
  downstream Met or a new stop cannot make a new protein.

check() applies this to the input isoforms before the update (a partial ORF is not complete, so the
isoform takes PASA's non-full-length path); check_models() to the transcripts the update added or
changed, whose CDS the update may have called again.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

from ..genome import Genome
from ..isoforms import Isoform, Orf
from ..models import AnnotationSet, Transcript
from . import checks as C
from .sites import POLYA, TSS, SiteIndex, five_prime, three_prime
from .structure import introns_of

log = logging.getLogger("pasa_xg")


def _order(t: Transcript):
    return t.exons if t.strand == "+" else list(reversed(t.exons))


def g2t(t: Transcript, pos: int) -> Optional[int]:
    """0-based transcript offset of genomic position pos, or None when pos is not exonic."""
    off = 0
    for a, b in _order(t):
        if a <= pos <= b:
            return off + (pos - a if t.strand == "+" else b - pos)
        off += b - a + 1
    return None


def cds_start(t: Transcript) -> int:
    return t.cds[0][0] if t.strand == "+" else t.cds[-1][1]


def cds_stop(t: Transcript) -> int:
    return t.cds[-1][1] if t.strand == "+" else t.cds[0][0]


def related_references(t: Transcript, ann: AnnotationSet) -> List[Tuple[Tuple[int, int, int], Transcript]]:
    """Same-strand coding reference transcripts sharing an intron with t (CDS overlap for an unspliced t),
    with a sort key: shared introns, CDS overlap, CDS length."""
    ints = set(introns_of(t.exons))
    out = []
    for g in ann.overlapping(t.contig, t.start, t.end, strand=t.strand):
        for r in g.transcripts.values():
            if not r.cds:
                continue
            shared = len(ints & set(introns_of(r.exons)))
            cds_ov = sum(max(0, min(b, d) - max(a, c) + 1) for a, b in t.exons for c, d in r.cds)
            if (ints and shared == 0) or cds_ov == 0:
                continue
            out.append(((shared, cds_ov, sum(b - a + 1 for a, b in r.cds)), r))
    return sorted(out, key=lambda kr: kr[0], reverse=True)


def _cds_phase(r: Transcript, pos: int) -> Optional[int]:
    """Position of genomic pos within its codon of r's CDS (0, 1, 2), or None outside the CDS."""
    off = 0
    for a, b in (r.cds if r.strand == "+" else list(reversed(r.cds))):
        if a <= pos <= b:
            off += pos - a if r.strand == "+" else b - pos
            return (off - r.cds_offset) % 3
        off += b - a + 1
    return None


def _positions(t: Transcript) -> Iterable[int]:
    for a, b in _order(t):
        yield from (range(a, b + 1) if t.strand == "+" else range(b, a - 1, -1))


def _read(seq: str, start: int, table: Dict[str, str]) -> Tuple[int, bool]:
    j = start
    while j + 3 <= len(seq):
        if table.get(seq[j:j + 3]) == "*":
            return j + 3, True
        j += 3
    return j, False


@dataclass
class Call:
    why: str                      # reference_start | reference_start_uorf | reference_frame_partial
    start: int = 0                # 0-based half-open ORF on the transcript sequence
    end: int = 0
    has_start: bool = True
    has_stop: bool = True


def decide(t: Transcript, refs, seq: str, table: Dict[str, str], p: C.Params, tss: bool, polya: bool,
           min_aa: int) -> Optional[Call]:
    """What the reference frame says about t's current CDS (None: keep it)."""
    if not t.cds or not refs:
        return None
    own = g2t(t, cds_start(t))
    if own is None:
        return None
    # the best-matching reference's start codon (any reference's 5'-most ATG fires on 2.5% of GENCODE's
    # own protein-coding transcripts, the best reference's on 1.0%)
    r_atg = g2t(t, cds_start(refs[0][1]))
    if r_atg is not None and seq[r_atg:r_atg + 3] == "ATG" and r_atg < own:
        s = r_atg
        e, has_stop = _read(seq, s, table)
        probe = t.copy(cds=t.segments_for(s + 1, e))
        if has_stop and (e - s) // 3 - 1 >= min_aa and C.nmd_distance(probe) <= p.NMD_STOP_TO_JUNCTION:
            return Call("reference_start", s, e, True, True)
        return Call("reference_start_uorf")
    r = refs[0][1]
    rs, re_ = cds_start(r), cds_stop(r)
    starts_inside = (t.start > rs) if t.strand == "+" else (t.end < rs)
    ends_inside = (t.end < re_) if t.strand == "+" else (t.start > re_)
    if not ((starts_inside and not tss) or (ends_inside and not polya)):
        return None
    s0 = g2t(t, rs)
    if s0 is not None and seq[s0:s0 + 3] == "ATG":
        start, has_start = s0, True
    else:
        start = next((i + (3 - ph) % 3 for i, pos in enumerate(_positions(t))
                      for ph in [_cds_phase(r, pos)] if ph is not None), None)
        if start is None:
            return None
        has_start = False
    end, has_stop = _read(seq, start, table)
    if has_start and has_stop:
        return None     # a complete ORF in the reference frame: nothing partial about it
    return Call("reference_frame_partial", start, end, has_start, has_stop)


def _apply(t: Transcript, c: Call) -> None:
    if c.why == "reference_start_uorf":
        t.cds = []
        t.attrs = {**t.attrs, "pasa_xg_orf": c.why}
        return
    t.cds = t.segments_for(c.start + 1, c.end)
    t.cds_offset = 0
    t.attrs = {**{k: v for k, v in t.attrs.items() if k not in ("cds_start_NF", "cds_end_NF")}, "pasa_xg_orf": c.why,
               **({"cds_start_NF": "true"} if not c.has_start else {}),
               **({"cds_end_NF": "true"} if not c.has_stop else {})}


def _label(c: Call) -> str:
    if c.why == "reference_start_uorf":
        return "noncoding: reference start opens a short / NMD ORF"
    if c.why == "reference_start":
        return "ORF from the reference start codon"
    return "partial in the reference frame"


def check(isos: Dict[str, Isoform], ann: AnnotationSet, genome: Genome, table: Dict[str, str],
          sites: Optional[SiteIndex], p: C.Params, min_aa: int = 100) -> Counter:
    """The check on input isoforms with an ORF (in place); counts by outcome."""
    from .noncoding import supported_ends   # avoids an import cycle
    n: Counter = Counter()
    for iso in isos.values():
        t = iso.tx
        if iso.orf is None or not t.cds:
            continue
        tss, polya = supported_ends(iso, sites)
        c = decide(t, related_references(t, ann), genome.transcript_seq(t), table, p, tss, polya, min_aa)
        if c is None:
            continue
        _apply(t, c)
        if c.why == "reference_start_uorf":
            iso.orf = None
        else:
            kind = ("complete" if c.has_stop else "3prime_partial") if c.has_start else \
                ("5prime_partial" if c.has_stop else "internal")
            iso.orf = Orf(f"{iso.id}.refframe", kind, c.start + 1, c.end, (c.end - c.start) // 3 - c.has_stop)
        n[_label(c)] += 1
    log.info("reference-frame check, isoforms: %s", dict(n))
    return n


def check_models(new: AnnotationSet, start: AnnotationSet, genome: Genome, table: Dict[str, str],
                 sites: Optional[SiteIndex], p: C.Params, min_aa: int = 100) -> Counter:
    """The check on the update's added transcripts and those whose CDS it changed (in place), against
    the starting annotation; a transcript left without CDS is labeled pasa_xg_coding 'none'."""
    old = {t.id: t for t in start.transcripts()}
    n: Counter = Counter()
    for t in new.transcripts():
        o = old.get(t.id)
        if not t.cds or (o is not None and o.cds == t.cds):
            continue    # reference models keep their CDS (UTR updates included)
        tss = sites is not None and sites.near(t.contig, t.strand, TSS, five_prime(t)) is not None
        polya = sites is not None and sites.near(t.contig, t.strand, POLYA, three_prime(t)) is not None
        c = decide(t, related_references(t, start), genome.transcript_seq(t), table, p, tss, polya, min_aa)
        if c is None:
            continue
        if o is not None:
            # a reference transcript keeps its own CDS (and its structure, when the CDS no longer fits)
            if not all(any(a <= x and y <= b for a, b in t.exons) for x, y in o.cds):
                t.exons = list(o.exons)
            t.cds, t.cds_offset = list(o.cds), o.cds_offset
            t.attrs = {**t.attrs, "pasa_xg_orf": "reference_cds_kept"}
            n["reference transcript: its CDS kept"] += 1
            continue
        _apply(t, c)
        if not t.cds:
            t.attrs = {**t.attrs, "pasa_xg_coding": "none"}
        n[_label(c)] += 1
    log.info("reference-frame check, updated models: %s", dict(n))
    return n


def keep_reference_cds(new: AnnotationSet, start: AnnotationSet, genome: Genome, table: Dict[str, str]) -> int:
    """A reference transcript whose protein the update left unchanged keeps its own CDS coordinates (an
    update can trim a partial model's incomplete codon or reset its start offset without changing the
    protein); only when the original CDS still lies in the transcript's exons."""
    old = {t.id: t for t in start.transcripts()}
    n = 0
    for t in new.transcripts():
        o = old.get(t.id)
        if o is None or not o.cds or not t.cds or (o.cds == t.cds and o.cds_offset == t.cds_offset):
            continue
        if not all(any(a <= x and y <= b for a, b in t.exons) for x, y in o.cds):
            continue
        if genome.protein(t, table).rstrip("*") == genome.protein(o, table).rstrip("*"):
            t.cds, t.cds_offset = list(o.cds), o.cds_offset
            n += 1
    log.info("reference transcripts with an unchanged protein given back their CDS: %d", n)
    return n


def sync_partial_tags(new: AnnotationSet, start: AnnotationSet, genome: Genome, table: Dict[str, str]) -> Counter:
    """cds_start_NF / cds_end_NF on every added or CDS-changed transcript from its own sequence (start
    codon, stop codon), in the attributes and in GENCODE's tag list, so a reference's stale NF tag goes
    when its CDS is completed and a partial CDS is always marked."""
    from .checks import tags
    old = {t.id: t for t in start.transcripts()}
    n: Counter = Counter()
    for t in new.transcripts():
        o = old.get(t.id)
        if o is not None and o.cds == t.cds:
            continue
        tl = [x for x in t.attrs.get("tag", "").split(",") if x and x not in ("cds_start_NF", "cds_end_NF")]
        attrs = {k: v for k, v in t.attrs.items() if k not in ("cds_start_NF", "cds_end_NF")}
        if t.cds:
            pep = genome.protein(t, table)
            same_start = o is not None and o.cds and cds_start(o) == cds_start(t) and "non_ATG_start" in tags(o)
            for k, ok in (("cds_start_NF", pep.startswith("M") or same_start), ("cds_end_NF", pep.endswith("*"))):
                if not ok:
                    attrs[k] = "true"
                    tl.append(k)
                    n[k] += 1
        if tl:
            attrs["tag"] = ",".join(tl)
        else:
            attrs.pop("tag", None)
        t.attrs = attrs
    log.info("partial CDS tags: %s", dict(n))
    return n
