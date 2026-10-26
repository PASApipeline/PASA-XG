"""Genome sequence access, reverse complement and translation."""

from __future__ import annotations

from typing import Dict, List, Tuple

import pysam

from .models import Transcript

_COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


_BASES = "TCAG"
_STANDARD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

# NCBI translation tables (amino acids in TCAG order), named as TransDecoder names
# them so one --genetic_code serves both. Start codons are taken to be ATG.
GENETIC_CODES: Dict[str, str] = {
    "Universal": _STANDARD,
    "Euplotid": "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "Ciliate": "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "Tetrahymena": "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "Candida": "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "Acetabularia": "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "Mitochondrial-Vertebrates": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
}


def codon_table(name: str = "Universal") -> Dict[str, str]:
    aas = GENETIC_CODES[name]
    table, i = {}, 0
    for a in _BASES:
        for b in _BASES:
            for c in _BASES:
                table[a + b + c] = aas[i]
                i += 1
    return table


def translate(seq: str, table: Dict[str, str]) -> str:
    seq = seq.upper()
    return "".join(table.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


class Genome:
    """Indexed FASTA; contigs are loaded once and cached."""

    def __init__(self, path: str) -> None:
        self.path = path
        self._fa = pysam.FastaFile(path)
        self._cache: Dict[str, str] = {}

    @property
    def contigs(self) -> List[str]:
        return list(self._fa.references)

    def length(self, contig: str) -> int:
        return self._fa.get_reference_length(contig)

    def contig(self, name: str) -> str:
        if name not in self._cache:
            self._cache[name] = self._fa.fetch(name).upper()
        return self._cache[name]

    def subseq(self, contig: str, start: int, end: int) -> str:
        return self.contig(contig)[start - 1:end]

    def spliced(self, contig: str, strand: str, spans: List[Tuple[int, int]]) -> str:
        seq = self.contig(contig)
        s = "".join(seq[a - 1:b] for a, b in sorted(spans))
        return revcomp(s) if strand == "-" else s

    def transcript_seq(self, t: Transcript) -> str:
        return self.spliced(t.contig, t.strand, t.exons)

    def cds_seq(self, t: Transcript) -> str:
        return self.spliced(t.contig, t.strand, t.cds)

    def protein(self, t: Transcript, table: Dict[str, str]) -> str:
        """Translation of the CDS from its first full codon (honors cds_offset)."""
        return translate(self.cds_seq(t)[t.cds_offset:], table) if t.cds else ""

    def splice_sites(self, t: Transcript) -> List[str]:
        """Donor..acceptor dinucleotides of each intron, on the transcript's strand."""
        seq = self.contig(t.contig)
        out = []
        for a, b in t.introns():
            d = seq[a - 1:a + 1] + ".." + seq[b - 2:b]
            if t.strand == "-":
                d = revcomp(seq[b - 2:b]) + ".." + revcomp(seq[a - 1:a + 1])
            out.append(d)
        return out
