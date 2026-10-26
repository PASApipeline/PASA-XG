"""Gene, transcript and annotation-set model.

Coordinates are 1-based and inclusive, as in GFF. Exons and CDS segments are
kept sorted by genomic position whatever the strand.
"""

from __future__ import annotations

import bisect
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

Span = Tuple[int, int]


def _merge_sorted(spans: Iterable[Span]) -> List[Span]:
    return sorted((int(a), int(b)) for a, b in spans)


@dataclass
class Transcript:
    id: str
    gene_id: str
    contig: str
    strand: str
    exons: List[Span]
    cds: List[Span] = field(default_factory=list)
    attrs: Dict[str, str] = field(default_factory=dict)
    # phase of the first CDS base in transcription order: a 5'-incomplete CDS
    # (e.g. GENCODE cds_start_NF) begins this many bases before its first full codon
    cds_offset: int = 0

    def __post_init__(self) -> None:
        if self.strand not in ("+", "-"):
            raise ValueError(f"{self.id}: strand must be + or -, got {self.strand!r}")
        self.exons = _merge_sorted(self.exons)
        self.cds = _merge_sorted(self.cds)

    @property
    def start(self) -> int:
        return self.exons[0][0]

    @property
    def end(self) -> int:
        return self.exons[-1][1]

    @property
    def is_coding(self) -> bool:
        return bool(self.cds)

    @property
    def num_exons(self) -> int:
        return len(self.exons)

    def introns(self) -> List[Span]:
        return [(self.exons[i][1] + 1, self.exons[i + 1][0] - 1) for i in range(len(self.exons) - 1)]

    def length(self) -> int:
        return sum(b - a + 1 for a, b in self.exons)

    def cds_length(self) -> int:
        return sum(b - a + 1 for a, b in self.cds)

    def cds_start(self) -> Optional[int]:
        """Genomic coordinate of the first CDS base in transcription order."""
        if not self.cds:
            return None
        return self.cds[0][0] if self.strand == "+" else self.cds[-1][1]

    def cds_end(self) -> Optional[int]:
        if not self.cds:
            return None
        return self.cds[-1][1] if self.strand == "+" else self.cds[0][0]

    def cds_phases(self) -> List[int]:
        """GFF phase of each CDS segment, in genomic order."""
        order = self.cds if self.strand == "+" else list(reversed(self.cds))
        phases, done = [], 0
        for a, b in order:
            phases.append((3 - (done - self.cds_offset) % 3) % 3)
            done += b - a + 1
        return phases if self.strand == "+" else list(reversed(phases))

    def utrs(self) -> Tuple[List[Span], List[Span]]:
        """(5' UTR, 3' UTR) segments in genomic order."""
        if not self.cds:
            return [], []
        lo, hi = self.cds[0][0], self.cds[-1][1]
        left, right = [], []
        for a, b in self.exons:
            if a < lo:
                left.append((a, min(b, lo - 1)))
            if b > hi:
                right.append((max(a, hi + 1), b))
        return (left, right) if self.strand == "+" else (right, left)

    def overlaps(self, other: "Transcript") -> bool:
        return (self.contig == other.contig and self.start <= other.end
                and other.start <= self.end)

    def exon_overlap(self, other: "Transcript") -> int:
        """Bases shared by the two transcripts' exons (same contig assumed)."""
        total, j = 0, 0
        for a, b in self.exons:
            while j < len(other.exons) and other.exons[j][1] < a:
                j += 1
            k = j
            while k < len(other.exons) and other.exons[k][0] <= b:
                total += max(0, min(b, other.exons[k][1]) - max(a, other.exons[k][0]) + 1)
                k += 1
        return total

    def transcript_to_genome(self, pos: int) -> int:
        """Genomic coordinate of 1-based transcript position `pos`."""
        order = self.exons if self.strand == "+" else list(reversed(self.exons))
        for a, b in order:
            n = b - a + 1
            if pos <= n:
                return a + pos - 1 if self.strand == "+" else b - pos + 1
            pos -= n
        raise ValueError(f"{self.id}: position beyond transcript length")

    def segments_for(self, t_start: int, t_end: int) -> List[Span]:
        """Genomic segments covering transcript positions t_start..t_end (1-based)."""
        out, offset = [], 0
        order = self.exons if self.strand == "+" else list(reversed(self.exons))
        for a, b in order:
            n = b - a + 1
            lo, hi = max(t_start, offset + 1), min(t_end, offset + n)
            if lo <= hi:
                if self.strand == "+":
                    out.append((a + lo - offset - 1, a + hi - offset - 1))
                else:
                    out.append((b - (hi - offset) + 1, b - (lo - offset) + 1))
            offset += n
        return sorted(out)

    def copy(self, **changes) -> "Transcript":
        d = dict(id=self.id, gene_id=self.gene_id, contig=self.contig, strand=self.strand,
                 exons=list(self.exons), cds=list(self.cds), attrs=dict(self.attrs),
                 cds_offset=self.cds_offset)
        d.update(changes)
        return Transcript(**d)


@dataclass
class Gene:
    id: str
    contig: str
    strand: str
    transcripts: Dict[str, Transcript] = field(default_factory=dict)
    attrs: Dict[str, str] = field(default_factory=dict)

    @property
    def start(self) -> int:
        return min(t.start for t in self.transcripts.values())

    @property
    def end(self) -> int:
        return max(t.end for t in self.transcripts.values())

    def add(self, t: Transcript) -> None:
        if t.contig != self.contig or t.strand != self.strand:
            raise ValueError(f"transcript {t.id} is not on gene {self.id}'s contig/strand")
        t.gene_id = self.id
        self.transcripts[t.id] = t

    def copy(self) -> "Gene":
        g = Gene(self.id, self.contig, self.strand, attrs=dict(self.attrs))
        for t in self.transcripts.values():
            g.transcripts[t.id] = t.copy()
        return g


class AnnotationSet:
    """Genes by id, with an overlap index per contig."""

    def __init__(self, genes: Iterable[Gene] = (), label: str = "input") -> None:
        self.label = label
        self.genes: Dict[str, Gene] = {}
        for g in genes:
            self.genes[g.id] = g
        self._index: Optional[Dict[str, Tuple[List[int], List[Tuple[int, int, str]]]]] = None

    def __len__(self) -> int:
        return len(self.genes)

    def __iter__(self) -> Iterator[Gene]:
        return iter(self.genes.values())

    def transcripts(self) -> Iterator[Transcript]:
        for g in self.genes.values():
            yield from g.transcripts.values()

    def add_gene(self, g: Gene) -> None:
        self.genes[g.id] = g
        self._index = None

    def remove_gene(self, gene_id: str) -> None:
        del self.genes[gene_id]
        self._index = None

    def copy(self, label: Optional[str] = None) -> "AnnotationSet":
        return AnnotationSet((g.copy() for g in self.genes.values()), label or self.label)

    def _build_index(self) -> None:
        by_contig = defaultdict(list)
        for g in self.genes.values():
            by_contig[g.contig].append((g.start, g.end, g.id))
        self._index = {}
        for contig, rows in by_contig.items():
            rows.sort()
            self._index[contig] = ([r[0] for r in rows], rows)
        # longest gene per contig bounds the backward scan
        self._max_len = {c: max(r[1] - r[0] + 1 for r in rows) for c, (_, rows) in self._index.items()}

    def overlapping(self, contig: str, start: int, end: int,
                    strand: Optional[str] = None) -> List[Gene]:
        if self._index is None:
            self._build_index()
        if contig not in self._index:
            return []
        starts, rows = self._index[contig]
        lo = bisect.bisect_left(starts, start - self._max_len[contig])
        hi = bisect.bisect_right(starts, end)
        out = []
        for s, e, gid in rows[lo:hi]:
            if s <= end and e >= start:
                g = self.genes[gid]
                if strand is None or g.strand == strand:
                    out.append(g)
        return out
