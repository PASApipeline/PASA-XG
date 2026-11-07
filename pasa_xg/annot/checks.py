"""Validation checks used by the annotation update (PASA's Gene_validator,
compare_updated_proteins, check_UTRs and check_num_UTR_exons)."""

from __future__ import annotations

import re
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
    # PASA-XG: additive update for curated annotations. No existing model is replaced,
    # restitched or split; evidence that conflicts with a model becomes an alternative
    # isoform (or a new gene). UTR extensions and compatible CDS refinements still apply,
    # and genes still merge under PASA's guards (none of them has full-length support,
    # so readthrough between expressed neighbours can't fuse them).
    ADDITIVE: bool = False
    # PASA's --TRUST_FL_STATUS: a compatible full-length isoform that the model extends
    # past is still full-length (LRAA trims isoform ends), not demoted to an EST.
    TRUST_FL_STATUS: bool = False
    # PASA-XG: no update may produce a model whose stop codon lies more than
    # NMD_STOP_TO_JUNCTION bases upstream of its last exon-exon junction.
    REJECT_NMD_TARGETS: bool = False
    NMD_STOP_TO_JUNCTION: int = 50
    # PASA-XG: an isoform that differs from a model (or from another new isoform) only at
    # its ends becomes a new isoform only when its protein differs ("protein"); with TSS /
    # polyA sites, the ends where it differs must also be supported. "all": PASA's
    # behavior, every such isoform that passes the protein checks is added.
    ALT_END_ISOFORMS: str = "protein"
    MIN_SAMPLE_FRACTION: float = 0.01  # noncoding pass: sample minimum is max(2, fraction x n_samples_total)
    # PASA-XG: when an update would change a reference transcript's CDS in place (its structure still
    # holding the reference CDS): "extend" lets it extend the CDS (the reference protein stays, in frame,
    # within the new one), "never" keeps the reference CDS, "all" is PASA's behavior. A change not allowed
    # is added as a new isoform carrying the other coding region. --additive implies "never".
    REFERENCE_CDS_UPDATES: str = "extend"


CONSENSUS = {"GT..AG", "GC..AG"}


def partial_flags(t: Transcript, genome: Genome, table: Dict[str, str]) -> Tuple[bool, bool]:
    """(5' partial, 3' partial) of a coding model, read off its protein."""
    pep = genome.protein(t, table)
    return (not pep.startswith("M"), not pep.endswith("*"))


def nmd_distance(t: Transcript) -> int:
    """Bases from the stop codon's last base to the last exon-exon junction, along the
    transcript (<= 0 when the stop codon is in the last exon or the model is unspliced)."""
    if t.num_exons < 2 or not t.cds:
        return 0
    exons = t.exons if t.strand == "+" else t.exons[::-1]
    stop = t.cds[-1][1] if t.strand == "+" else t.cds[0][0]
    pos = 0
    for a, b in exons:
        if a <= stop <= b:
            pos += (stop - a + 1) if t.strand == "+" else (b - stop + 1)
            break
        pos += b - a + 1
    junction = t.length() - (exons[-1][1] - exons[-1][0] + 1)
    return junction - pos


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
        if p.REJECT_NMD_TARGETS and pep.endswith("*"):
            d = nmd_distance(t)
            if d > p.NMD_STOP_TO_JUNCTION:
                errs.append(f"NMD target: stop codon {d} bases upstream of the last splice junction")
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


def cds_extends(old: Transcript, new: Transcript, genome: Genome, table: Dict[str, str]) -> bool:
    """new's CDS extends old's: every base of old's CDS is in new's, and old's protein lies within
    new's (so the frame is kept). A shorter or reframed CDS is not an extension."""
    if not old.cds or not new.cds:
        return False
    ob = {x for a, b in old.cds for x in range(a, b + 1)}
    nb = {x for a, b in new.cds for x in range(a, b + 1)}
    if not ob <= nb:
        return False
    return genome.protein(old, table).rstrip("*") in genome.protein(new, table).rstrip("*")


IG_TR_TYPE = re.compile(r"^(IG|TR)_|^[VDJ]_segment|^C_region")   # GENCODE gene_type / RefSeq gene_biotype


def is_ig_tr(attrs: Dict[str, str]) -> bool:
    """Immunoglobulin / T-cell receptor gene segment (or segment pseudogene). Rearranged (V(D)J) BCR / TCR
    transcripts join segments somatically, so they are not germline isoforms of these genes."""
    return any(IG_TR_TYPE.match(attrs.get(k, "")) for k in ("gene_type", "gene_biotype"))


FROZEN_CDS_TAGS = {"seleno", "non_ATG_start"}   # GENCODE tags: the CDS can't be re-derived from the sequence


def tags(t: Transcript) -> set:
    return {x for x in t.attrs.get("tag", "").split(",") if x}


def compare_proteins(old: Transcript, new: Transcript, genome: Genome, table: Dict[str, str],
                     p: Params, fl_alt_splice: bool = False,
                     new_partial5: bool = False, new_partial3: bool = False) -> Tuple[bool, str]:
    """PASA's compare_updated_proteins: the new model must validate, keep most of
    the old protein's length, and align to it closely."""
    old_pep = genome.protein(old, table)
    new_pep = genome.protein(new, table)
    if not old_pep or not new_pep:
        return False, "no protein for comparison"
    if tags(old) & FROZEN_CDS_TAGS and new_pep != old_pep:
        # a selenoprotein's UGA or a non-ATG start: an ORF call on the sequence would change the protein
        return False, f"model's CDS is kept ({', '.join(sorted(tags(old) & FROZEN_CDS_TAGS))})"
    old_complete = (old_pep.startswith("M") or "non_ATG_start" in tags(old)) and old_pep.endswith("*")
    if old_complete and not (new_pep.startswith("M") and new_pep.endswith("*")):
        return False, "would replace a complete protein with a partial one"
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
