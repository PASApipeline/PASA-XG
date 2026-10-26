"""LRAA isoforms: structure plus expression and end support."""

from __future__ import annotations

import pickle
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import pandas as pd

from .gff import read_annotation
from .models import AnnotationSet, Transcript


@dataclass
class Orf:
    orf_id: str
    orf_type: str          # complete, 5prime_partial, 3prime_partial, internal
    t_start: int           # 1-based transcript coords of the CDS (incl. stop codon)
    t_end: int
    aa_len: int
    score: float = 0.0

    @property
    def complete(self) -> bool:
        return self.orf_type == "complete"


@dataclass
class Isoform:
    tx: Transcript
    uniq_reads: float = 0.0
    all_reads: float = 0.0
    isoform_fraction: float = 0.0
    tpm: float = 0.0
    tss: bool = False
    polya: bool = False
    orf: Optional[Orf] = None
    attrs: Dict[str, str] = field(default_factory=dict)

    @property
    def id(self) -> str:
        return self.tx.id

    @property
    def full_length(self) -> bool:
        """PASA's full-length path: a complete ORF."""
        return self.orf is not None and self.orf.complete


def _bool(v: Optional[str]) -> bool:
    return str(v).lower() == "true"


def load_lraa(gtf_path: str, expr_path: str) -> Dict[str, Isoform]:
    ann: AnnotationSet = read_annotation(gtf_path, label="lraa")
    expr = pd.read_csv(expr_path, sep="\t", comment="#")
    expr = expr.set_index("transcript_id")
    out: Dict[str, Isoform] = {}
    for t in ann.transcripts():
        iso = Isoform(t, tss=_bool(t.attrs.get("TSS")), polya=_bool(t.attrs.get("PolyA")))
        if t.id in expr.index:
            r = expr.loc[t.id]
            iso.uniq_reads = float(r["uniq_reads"])
            iso.all_reads = float(r["all_reads"])
            iso.isoform_fraction = float(r["isoform_fraction"])
            iso.tpm = float(r["TPM"])
        out[t.id] = iso
    return out


def table(isoforms: Dict[str, Isoform]) -> pd.DataFrame:
    rows = []
    for iso in isoforms.values():
        t = iso.tx
        rows.append({
            "transcript_id": t.id, "gene_id": t.gene_id, "contig": t.contig, "strand": t.strand,
            "start": t.start, "end": t.end, "num_exons": t.num_exons, "length": t.length(),
            "exons": ",".join(f"{a}-{b}" for a, b in t.exons),
            "uniq_reads": iso.uniq_reads, "all_reads": iso.all_reads,
            "isoform_fraction": iso.isoform_fraction, "TPM": iso.tpm,
            "TSS": iso.tss, "PolyA": iso.polya,
            "orf_type": iso.orf.orf_type if iso.orf else None,
            "orf_aa_len": iso.orf.aa_len if iso.orf else None,
        })
    return pd.DataFrame(rows)


def save(isoforms: Dict[str, Isoform], path: str) -> None:
    with open(path, "wb") as fh:
        pickle.dump(isoforms, fh)


def load(path: str) -> Dict[str, Isoform]:
    with open(path, "rb") as fh:
        return pickle.load(fh)
