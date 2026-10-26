"""Validation checks used by the annotation update (PASA's Gene_validator,
compare_updated_proteins, check_UTRs and check_num_UTR_exons)."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, List, Sequence, Tuple

from Bio import Align
from Bio.Align import substitution_matrices

from ..genome import Genome, translate
from ..models import Span, Transcript


@dataclass
class Params:
    """PASA's annotation-compare thresholds, same names and defaults."""
    MIN_PERCENT_OVERLAP: float = 50
    MIN_PERCENT_PROT_CODING: float = 40
    MIN_PERID_PROT_COMPARE: float = 70
    MIN_PERCENT_LENGTH_FL_COMPARE: float = 70
    MIN_PERCENT_LENGTH_NONFL_COMPARE: float = 70
    MIN_FL_ORF_SIZE: int = 0          # aa; when set, overrides MIN_PERCENT_PROT_CODING
    MIN_PERCENT_ALIGN_LENGTH: float = 70
    MIN_PERCENT_OVERLAP_GENE_REPLACE: float = 80
    STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE: bool = False
    MAX_UTR_EXONS: int = 2
    ALLOW_ATAC: bool = True
    # PASA-XG: a model whose own structure the isoforms support is never replaced
    # (PASA status 8) or restitched (16); the new structure becomes an alternative
    # isoform instead. False restores PASA's behavior.
    PROTECT_SUPPORTED_MODELS: bool = True


CONSENSUS = {"GT..AG", "GC..AG"}


def partial_flags(t: Transcript, genome: Genome, table: Dict[str, str]) -> Tuple[bool, bool]:
    """(5' partial, 3' partial) of a coding model, read off its protein."""
    pep = genome.protein(t, table)
    return (not pep.startswith("M"), not pep.endswith("*"))


def validate_gene(t: Transcript, genome: Genome, table: Dict[str, str], p: Params,
                  partial5: bool = False, partial3: bool = False) -> List[str]:
    """Structural and coding checks; an empty list means the model validates."""
    errs: List[str] = []
    if not t.cds:
        return ["no CDS"]
    if not all(any(ea <= a and b <= eb for ea, eb in t.exons) for a, b in t.cds):
        errs.append("CDS outside exons")
    ok_sites = CONSENSUS | ({"AT..AC"} if p.ALLOW_ATAC else set())
    bad = [s for s in genome.splice_sites(t) if s not in ok_sites]
    if bad:
        errs.append("non-consensus splice sites " + ",".join(bad))
    if t.start < 1 or t.end > genome.length(t.contig):
        errs.append("extends beyond the contig")
    pep = genome.protein(t, table)
    if pep.count("*") > 1:
        errs.append(f"{pep.count('*')} stops in protein")
    else:
        if not partial5 and not pep.startswith("M"):
            errs.append("protein does not start with M")
        if not partial3 and not pep.endswith("*"):
            errs.append("protein lacks a stop codon")
    return errs


@lru_cache(maxsize=1)
def _aligner() -> Align.PairwiseAligner:
    a = Align.PairwiseAligner()
    a.mode = "local"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -10
    a.extend_gap_score = -1
    return a


def align_proteins(a: str, b: str) -> Tuple[float, int]:
    """(percent identity, alignment length in residues) of the best local alignment."""
    a, b = a.rstrip("*").replace("*", "X"), b.rstrip("*").replace("*", "X")
    if not a or not b:
        return 0.0, 0
    aln = _aligner().align(a, b)[0]
    c = aln.counts()
    length = c.aligned + c.internal_gaps
    return (100.0 * c.identities / length if length else 0.0), length


def compare_proteins(old: Transcript, new: Transcript, genome: Genome, table: Dict[str, str],
                     p: Params, fl_alt_splice: bool = False,
                     new_partial5: bool = False, new_partial3: bool = False) -> Tuple[bool, str]:
    """PASA's compare_updated_proteins: the new model must validate, keep most of
    the old protein's length, and align to it closely."""
    old_pep = genome.protein(old, table)
    new_pep = genome.protein(new, table)
    if not old_pep or not new_pep:
        return False, "no protein for comparison"
    errs = validate_gene(new, genome, table, p, new_partial5, new_partial3)
    if errs:
        return False, "; ".join(errs)
    pct_len = 100.0 * len(new_pep) / len(old_pep)
    min_len = p.MIN_PERCENT_LENGTH_FL_COMPARE if fl_alt_splice else p.MIN_PERCENT_LENGTH_NONFL_COMPARE
    if pct_len < min_len:
        return False, f"new protein is {pct_len:.0f}% of the original length (< {min_len:g}%)"
    per_id, overlap = align_proteins(old_pep, new_pep)
    shorter = min(len(old_pep), len(new_pep)) - (1 if old_pep.endswith("*") and new_pep.endswith("*") else 0)
    pct_overlap = 100.0 * overlap / max(1, shorter)
    if per_id < p.MIN_PERID_PROT_COMPARE or pct_overlap < p.MIN_PERCENT_ALIGN_LENGTH:
        return False, (f"new protein aligns at {per_id:.0f}% identity over {pct_overlap:.0f}% "
                       "of the shorter protein")
    return True, f"protein aligns at {per_id:.0f}% identity over {min(100, pct_overlap):.0f}%"


def utr_exon_count(t: Transcript) -> int:
    """Largest run of adjacent exons without CDS."""
    best = run = 0
    for a, b in t.exons:
        if any(ca <= b and cb >= a for ca, cb in t.cds):
            run = 0
        else:
            run += 1
            best = max(best, run)
    return best


def check_utr_exons(t: Transcript, p: Params) -> Tuple[bool, str]:
    n = utr_exon_count(t)
    if n > p.MAX_UTR_EXONS:
        return False, f"{n} adjacent UTR exons (> {p.MAX_UTR_EXONS})"
    return True, ""


def check_false_utrs(evidence_span: Span, old: Transcript, new: Transcript) -> Tuple[bool, str]:
    """PASA's check_UTRs: after stitching, a CDS shift that leaves the new model with
    a UTR the original lacked, extending past the evidence, is an artifact."""
    if not old.cds or not new.cds:
        return True, ""
    om_l, om_r = old.cds[0][0], old.cds[-1][1]
    nm_l, nm_r = new.cds[0][0], new.cds[-1][1]
    if (om_l, om_r) == (nm_l, nm_r):
        return True, ""
    five, three = old.utrs()
    left_utr, right_utr = (five, three) if old.strand == "+" else (three, five)
    al, ar = evidence_span
    msgs = []
    if not left_utr and om_l != nm_l and om_l == new.start and om_l < al:
        msgs.append("left UTR not supported by the transcript")
    if not right_utr and om_r != nm_r and om_r == new.end and om_r > ar:
        msgs.append("right UTR not supported by the transcript")
    return (not msgs), "; ".join(msgs)


def span_overlap(a: Span, b: Span) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


def cds_span(t: Transcript) -> Span:
    return (t.cds[0][0], t.cds[-1][1]) if t.cds else (t.start, t.end)


def same_exons(a: Sequence[Span], b: Sequence[Span]) -> bool:
    return list(a) == list(b)
