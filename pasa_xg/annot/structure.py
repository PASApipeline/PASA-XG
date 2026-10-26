"""Transcript-structure operations behind the annotation update.

Ports of PASA's pairwise compatibility test (PASA_alignment_assembler and
Gene_obj_alignment_assembler for two members), CDNA_stitcher, and the longest-ORF
gene construction of Exons_to_geneobj.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from ..genome import Genome, translate
from ..models import Span, Transcript

FUZZ = 20  # PASA's fuzzlength: a terminal exon may run this far into an intron


def introns_of(exons: Sequence[Span]) -> List[Span]:
    return [(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)]


def exons_from(start: int, end: int, introns: Sequence[Span]) -> List[Span]:
    out, pos = [], start
    for a, b in sorted(introns):
        out.append((pos, a - 1))
        pos = b + 1
    out.append((pos, end))
    return out


def _trim_ends(x: Sequence[Span], y_introns: Sequence[Span], fuzz: int) -> Optional[Tuple[int, int]]:
    """x's span, with a terminal exon that runs at most `fuzz` bases into one of y's
    introns trimmed back to the splice site; None if x conflicts with y's introns."""
    start, end = x[0][0], x[-1][1]
    x_introns = set(introns_of(x))
    for a, b in y_introns:
        if b < start or a > end or (a, b) in x_introns:
            continue
        if a <= start <= b:  # x starts inside y's intron
            if x[0][1] <= b or b - start + 1 > fuzz:
                return None
            start = b + 1
        elif a <= end <= b:  # x ends inside y's intron
            if x[-1][0] >= a or end - a + 1 > fuzz:
                return None
            end = a - 1
        else:  # inside x's span, and x lacks it: an exon or another intron conflicts
            return None
    return start, end


def merge_compatible(x: Sequence[Span], y: Sequence[Span], fuzz: int = FUZZ) -> Optional[List[Span]]:
    """The single structure two overlapping transcripts assemble into, or None if
    they are incompatible (PASA: they would not assemble into one)."""
    x, y = sorted(x), sorted(y)
    if x[-1][1] < y[0][0] or y[-1][1] < x[0][0]:
        return None
    if not _exons_overlap(x, y):
        return None
    tx = _trim_ends(x, introns_of(y), fuzz)
    ty = _trim_ends(y, introns_of(x), fuzz)
    if tx is None or ty is None:
        return None
    start, end = min(tx[0], ty[0]), max(tx[1], ty[1])
    introns = sorted(set(introns_of(x)) | set(introns_of(y)))
    for i in range(1, len(introns)):
        if introns[i][0] <= introns[i - 1][1] + 1:
            return None  # overlapping or abutting distinct introns
    introns = [iv for iv in introns if iv[0] > start and iv[1] < end]
    return exons_from(start, end, introns)


def _exons_overlap(x: Sequence[Span], y: Sequence[Span]) -> bool:
    j = 0
    for a, b in x:
        while j < len(y) and y[j][1] < a:
            j += 1
        if j < len(y) and y[j][0] <= b:
            return True
    return False


def compatible(x: Sequence[Span], y: Sequence[Span], fuzz: int = FUZZ) -> bool:
    return merge_compatible(x, y, fuzz) is not None


def stitch(template: Sequence[Span], thread: Sequence[Span], fuzz: int = FUZZ) -> List[Span]:
    """PASA's CDNA_stitcher: replace the template's structure across the thread's
    span with the thread's structure.

    Segments carry splice-junction flags: an alignment's first exon has no left
    junction and its last no right junction. The thread's terminal exons merge with
    the template exons they overlap; a thread end within `fuzz` of a template splice
    site snaps to it, keeping the template's exons beyond.
    """
    T = [[a, b, i > 0, i < len(template) - 1] for i, (a, b) in enumerate(sorted(template))]
    H = [[a, b, i > 0, i < len(thread) - 1] for i, (a, b) in enumerate(sorted(thread))]
    single = len(H) == 1
    out: List[List] = []

    def anchor(lend: int, rend: int, reverse: bool) -> int:
        order = range(len(T) - 1, -1, -1) if reverse else range(len(T))
        for i in order:
            if T[i][0] < rend and T[i][1] > lend:
                return i
        return -1

    # first threaded segment
    f = H[0]
    k = anchor(f[0], f[1], False)
    if k > -1:
        al, ar, aleft, aright = T[k]
        new_l, new_r = min(al, f[0]), max(ar, f[1])
        has_left = aleft
        if aleft and al - f[0] <= fuzz:
            new_l = al
        elif f[0] < al:
            has_left = f[2]
        has_right = f[3]
        if f[3]:
            new_r = f[1]
        elif aright and f[1] - ar <= fuzz:
            new_r = ar
            has_right = aright
        if has_left:
            out.extend([list(s) for s in T[:k]])
        out.append([new_l, new_r, has_left, has_right])
    else:
        out.append(list(f))
    # internal threaded segments
    out.extend([list(s) for s in H[1:-1]])
    # last threaded segment
    last = H[-1]
    k = anchor(last[0], last[1], True)
    if k > -1:
        al, ar, aleft, aright = T[k]
        has_left, has_right = last[2], aright
        new_l, new_r = min(al, last[0]), max(ar, last[1])
        if last[2]:
            new_l = last[0]
        elif aleft and al - last[0] <= fuzz:
            new_l = al
            has_left = aleft
        if aright and last[1] - ar <= fuzz:
            new_r = ar
        else:
            has_right = last[3]
        if single:
            out[-1][1] = new_r
            out[-1][3] = has_right
        else:
            out.append([new_l, new_r, has_left, has_right])
    elif not single:
        out.append(list(last))
    # template exons beyond the thread, if the stitched structure splices onward
    if out[-1][3]:
        l, r = out[-1][0], out[-1][1]
        k = -1
        for i in range(len(T) - 1, -1, -1):
            if T[i][1] > l and T[i][0] < r:
                k = i
                break
        if k > -1:
            out.extend([list(s) for s in T[k + 1:]])
    exons = sorted((a, b) for a, b, _, _ in out)
    merged: List[List[int]] = []
    for a, b in exons:  # guard against touching segments
        if merged and a <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return [(a, b) for a, b in merged]


def longest_orf(seq: str, table: Dict[str, str], allow_5partial: bool = False,
                allow_3partial: bool = False) -> Optional[Tuple[int, int]]:
    """Longest forward-strand ORF, 1-based (start, end) including the stop codon.

    Complete ORFs run ATG..stop; a 5' partial ORF may start at the first codon of a
    frame, a 3' partial ORF may run off the end (PASA's Longest_orf)."""
    seq = seq.upper()
    best: Optional[Tuple[int, int]] = None

    def keep(s: int, e: int) -> None:
        nonlocal best
        if best is None or e - s > best[1] - best[0]:
            best = (s, e)

    for frame in range(3):
        start_idx: Optional[int] = 0 if allow_5partial else None
        ncodons = (len(seq) - frame) // 3
        for ci in range(ncodons):
            codon = seq[frame + 3 * ci: frame + 3 * ci + 3]
            if start_idx is None and codon == "ATG":
                start_idx = ci
            if table.get(codon) == "*":
                if start_idx is not None:
                    keep(frame + 3 * start_idx + 1, frame + 3 * ci + 3)
                start_idx = None
        if start_idx is not None and allow_3partial and ncodons > start_idx:
            keep(frame + 3 * start_idx + 1, frame + 3 * ncodons)
    return best


def gene_from_exons(exons: Sequence[Span], strand: str, contig: str, genome: Genome,
                    table: Dict[str, str], tid: str = "tmp", gene_id: str = "tmp",
                    allow_5partial: bool = False, allow_3partial: bool = False) -> Transcript:
    """A transcript over `exons` with its longest ORF as the CDS."""
    t = Transcript(tid, gene_id, contig, strand, list(exons))
    orf = longest_orf(genome.transcript_seq(t), table, allow_5partial, allow_3partial)
    if orf:
        t.cds = t.segments_for(*orf)
    return t


def protein_of(t: Transcript, genome: Genome, table: Dict[str, str]) -> str:
    return translate(genome.cds_seq(t), table) if t.cds else ""
