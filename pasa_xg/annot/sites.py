"""TSS and polyA sites for the annotation update (docs/TSS_POLYA_DESIGN.md).

Two sources: LRAA's own calls on isoform ends (the TSS / PolyA flags), and known sites
given as BED files (CAGE peaks, PolyASite clusters, ...). An isoform end is supported when
LRAA called it or a known site lies within the window; a model end is anchored when a known
site or a supported isoform end lies within the window.
"""

from __future__ import annotations

import bisect
import gzip
import os
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Tuple

from ..models import Span, Transcript

TSS, POLYA = "tss", "polya"


def bare(contig: str) -> str:
    """Contig names compared without a 'chr' prefix (PolyASite writes '22', GENCODE 'chr22')."""
    return contig[3:] if contig.startswith("chr") else contig


def kind_at(strand: str, left: bool) -> str:
    """Which site a transcript's left or right end is."""
    return TSS if left == (strand == "+") else POLYA


def end_of(exons: List[Span], left: bool) -> int:
    return exons[0][0] if left else exons[-1][1]


def five_prime(t: Transcript) -> int:
    return t.start if t.strand == "+" else t.end


def three_prime(t: Transcript) -> int:
    return t.end if t.strand == "+" else t.start


class SiteIndex:
    """Sites per (contig, strand, kind), as intervals with a source label."""

    def __init__(self, window: int = 50) -> None:
        self.window = window
        self._raw: Dict[Tuple[str, str, str], List[Tuple[int, int, str]]] = defaultdict(list)
        self._starts: Dict[Tuple[str, str, str], List[int]] = {}
        self._span: Dict[Tuple[str, str, str], int] = {}

    def __len__(self) -> int:
        return sum(len(v) for v in self._raw.values())

    def add(self, contig: str, strand: str, kind: str, start: int, end: int, source: str) -> None:
        self._raw[(bare(contig), strand, kind)].append((start, end, source))
        self._starts.clear()

    def _ready(self) -> None:
        if self._starts or not self._raw:
            return
        for k, v in self._raw.items():
            v.sort()
            self._starts[k] = [a for a, _, _ in v]
            self._span[k] = max(b - a for a, b, _ in v)

    def near(self, contig: str, strand: str, kind: str, pos: int,
             sources: Optional[Iterable[str]] = None) -> Optional[str]:
        """Source of the closest site within the window of pos, or None. `sources`
        restricts to labels starting with one of the given prefixes."""
        self._ready()
        k = (bare(contig), strand, kind)
        v = self._raw.get(k)
        if not v:
            return None
        w, starts = self.window, self._starts[k]
        i = bisect.bisect_right(starts, pos + w)
        best, best_d = None, w + 1
        lo = pos - w - self._span[k]
        while i > 0 and starts[i - 1] >= lo:
            i -= 1
            a, b, src = v[i]
            if sources is not None and not any(src.startswith(s) for s in sources):
                continue
            d = 0 if a <= pos <= b else min(abs(pos - a), abs(pos - b))
            if d < best_d:
                best, best_d = src, d
        return best


def load_bed(index: SiteIndex, path: str, kind: str) -> int:
    """BED6 (0-based start, strand in column 6). Sites without a strand are skipped."""
    name = os.path.basename(path)
    for suffix in (".gz", ".bed"):
        name = name[:-len(suffix)] if name.endswith(suffix) else name
    n = 0
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[5] not in ("+", "-"):
                continue
            index.add(f[0], f[5], kind, int(f[1]) + 1, int(f[2]), "known:" + name)
            n += 1
    return n


def add_isoform_ends(index: SiteIndex, isoforms) -> int:
    """LRAA-called ends of the isoforms (Isoform.tss / .polya)."""
    n = 0
    for i in isoforms:
        t = i.tx
        if i.tss:
            index.add(t.contig, t.strand, TSS, five_prime(t), five_prime(t), "lraa")
            n += 1
        if i.polya:
            index.add(t.contig, t.strand, POLYA, three_prime(t), three_prime(t), "lraa")
            n += 1
    return n


def iso_support(index: Optional[SiteIndex], iso, kind: str) -> Optional[str]:
    """Why an isoform end counts as supported ('lraa', 'known:<bed>'), or None."""
    if kind == TSS and iso.tss:
        return "lraa"
    if kind == POLYA and iso.polya:
        return "lraa"
    if index is None:
        return None
    t = iso.tx
    pos = five_prime(t) if kind == TSS else three_prime(t)
    return index.near(t.contig, t.strand, kind, pos, sources=("known:",))
