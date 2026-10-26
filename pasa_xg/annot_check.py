"""Checks on reference gene models: is each CDS a sound, translatable ORF?"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

from .genome import Genome, codon_table, translate
from .models import AnnotationSet, Transcript

CANONICAL = {"GT..AG", "GC..AG", "AT..AC"}


@dataclass
class CdsCheck:
    transcript_id: str
    gene_id: str
    coding: bool
    cds_in_exons: bool = True
    frame_ok: bool = True
    has_start: bool = True
    has_stop: bool = True
    internal_stops: int = 0
    noncanonical_introns: int = 0

    @property
    def ok(self) -> bool:
        return (not self.coding) or (self.cds_in_exons and self.frame_ok and self.internal_stops == 0)

    def row(self) -> Dict[str, object]:
        d = dict(self.__dict__)
        d["ok"] = self.ok
        return d


def protein(t: Transcript, genome: Genome, table: Dict[str, str]) -> str:
    return genome.protein(t, table)


def check_transcript(t: Transcript, genome: Genome, table: Dict[str, str]) -> CdsCheck:
    c = CdsCheck(t.id, t.gene_id, t.is_coding)
    c.noncanonical_introns = sum(1 for s in genome.splice_sites(t) if s not in CANONICAL)
    if not t.is_coding:
        return c
    c.cds_in_exons = all(any(ea <= a and b <= eb for ea, eb in t.exons) for a, b in t.cds)
    pep = protein(t, genome, table)
    c.has_start = pep.startswith("M")
    c.has_stop = pep.endswith("*")
    # a 3'-incomplete CDS (no stop) need not end on a codon boundary
    c.frame_ok = (not c.has_stop) or (t.cds_length() - t.cds_offset) % 3 == 0
    c.internal_stops = pep.rstrip("*").count("*")
    return c


def check_annotation(ann: AnnotationSet, genome: Genome, code: str = "Universal") -> List[CdsCheck]:
    table = codon_table(code)
    return [check_transcript(t, genome, table) for t in ann.transcripts()]
