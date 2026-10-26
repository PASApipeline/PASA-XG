"""Read GTF/GFF3 gene annotations; write GFF3, GTF and BED12."""

from __future__ import annotations

import gzip
import re
from collections import OrderedDict, defaultdict
from typing import Dict, Iterable, List, Optional, TextIO, Tuple
from urllib.parse import quote, unquote

from .models import AnnotationSet, Gene, Transcript

TRANSCRIPT_TYPES = {"mRNA", "transcript", "ncRNA", "lnc_RNA", "lncRNA", "pseudogenic_transcript",
                    "primary_transcript", "miRNA", "snRNA", "snoRNA", "rRNA", "tRNA"}


def _open(path: str) -> TextIO:
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def sniff_format(path: str) -> str:
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##gff-version 3"):
                return "gff3"
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            return "gtf" if re.search(r'\w+ "[^"]*"', f[8]) else "gff3"
    raise ValueError(f"{path}: no features found")


def _gff3_attrs(s: str) -> Dict[str, str]:
    out = {}
    for part in s.strip().strip(";").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = unquote(v.strip())
    return out


def _gtf_attrs(s: str) -> Dict[str, str]:
    return {m.group(1): m.group(2) for m in re.finditer(r'(\S+) "([^"]*)"', s)}


def read_annotation(path: str, label: str = "input") -> AnnotationSet:
    fmt = sniff_format(path)
    return _read_gtf(path, label) if fmt == "gtf" else _read_gff3(path, label)


def _assemble(tx_rows, gene_of, gene_attrs, tx_attrs, label) -> AnnotationSet:
    genes: Dict[str, Gene] = OrderedDict()
    for tid, rows in tx_rows.items():
        contig, strand = rows["loc"]
        if not rows["exon"]:
            if not rows["cds"]:
                continue
            rows["exon"] = list(rows["cds"])  # CDS-only models: exons = CDS
        gid = gene_of.get(tid, tid)
        offset = 0
        if rows["phase"]:  # phase of the first CDS row in transcription order
            first = (min if strand == "+" else max)(rows["phase"], key=lambda sp: sp[0])
            offset = rows["phase"][first]
        t = Transcript(tid, gid, contig, strand, _merge_adjacent(rows["exon"]),
                       _merge_adjacent(rows["cds"]), tx_attrs.get(tid, {}), cds_offset=offset)
        g = genes.get(gid)
        if g is None:
            g = genes[gid] = Gene(gid, contig, strand, attrs=gene_attrs.get(gid, {}))
        if g.contig != contig or g.strand != strand:
            raise ValueError(f"gene {gid} has transcripts on different contigs/strands")
        g.add(t)
    return AnnotationSet(genes.values(), label)


def _merge_adjacent(spans: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Sort and join touching or overlapping spans (GTF stop_codon next to CDS,
    UTR and CDS rows describing one exon)."""
    out: List[List[int]] = []
    for a, b in sorted(spans):
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def _read_gtf(path: str, label: str) -> AnnotationSet:
    tx_rows = OrderedDict()
    gene_of, gene_attrs, tx_attrs = {}, {}, {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            a = _gtf_attrs(f[8])
            if ftype == "gene":
                gene_attrs[a.get("gene_id")] = {k: v for k, v in a.items() if k != "gene_id"}
                continue
            tid = a.get("transcript_id")
            if not tid:
                continue
            rows = tx_rows.setdefault(tid, {"loc": (f[0], f[6]), "exon": [], "cds": [], "phase": {}})
            gene_of[tid] = a.get("gene_id", tid)
            span = (int(f[3]), int(f[4]))
            if ftype == "exon":
                rows["exon"].append(span)
            elif ftype in ("CDS", "stop_codon"):  # GTF CDS excludes the stop codon
                rows["cds"].append(span)
                if ftype == "CDS" and f[7].isdigit():
                    rows["phase"][span] = int(f[7])
            elif ftype == "transcript":
                tx_attrs[tid] = {k: v for k, v in a.items() if k not in ("gene_id", "transcript_id")}
    return _assemble(tx_rows, gene_of, gene_attrs, tx_attrs, label)


def _read_gff3(path: str, label: str) -> AnnotationSet:
    parent_of: Dict[str, str] = {}
    types: Dict[str, str] = {}
    feats: Dict[str, Dict[str, str]] = {}
    locs: Dict[str, Tuple[str, str]] = {}
    children = defaultdict(list)  # parent id -> [(type, span)]
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = _gff3_attrs(f[8])
            ftype, span = f[2], (int(f[3]), int(f[4]))
            if ftype in ("exon", "CDS", "stop_codon"):
                for p in a.get("Parent", "").split(","):
                    if p:
                        children[p].append(("exon" if ftype == "exon" else "cds", span,
                                            int(f[7]) if ftype == "CDS" and f[7].isdigit() else None))
                continue
            fid = a.get("ID")
            if not fid:
                continue
            types[fid], feats[fid], locs[fid] = ftype, a, (f[0], f[6])
            if "Parent" in a:
                parent_of[fid] = a["Parent"].split(",")[0]
    tx_rows, gene_of, gene_attrs, tx_attrs = OrderedDict(), {}, {}, {}
    for fid, ftype in types.items():
        is_tx = ftype in TRANSCRIPT_TYPES or (fid in children and fid in parent_of)
        if not is_tx:
            continue
        rows = {"loc": locs[fid], "exon": [], "cds": [], "phase": {}}
        for kind, span, phase in children.get(fid, []):
            rows[kind].append(span)
            if phase is not None:
                rows["phase"][span] = phase
        tx_rows[fid] = rows
        gid = parent_of.get(fid, fid)
        gene_of[fid] = gid
        tx_attrs[fid] = {k: v for k, v in feats[fid].items() if k not in ("ID", "Parent")}
        if gid in feats:
            gene_attrs[gid] = {k: v for k, v in feats[gid].items() if k not in ("ID",)}
    return _assemble(tx_rows, gene_of, gene_attrs, tx_attrs, label)


# ---------------------------------------------------------------- writers

def _fmt_gff3_attrs(d: Dict[str, str]) -> str:
    return ";".join(f"{k}={quote(str(v), safe=' :|.,-_/()+')}" for k, v in d.items())


def write_gff3(genes: Iterable[Gene], fh: TextIO, source: str = "PASA-XG") -> None:
    fh.write("##gff-version 3\n")
    for g in sorted(genes, key=lambda g: (g.contig, g.start, g.id)):
        fh.write(f"{g.contig}\t{source}\tgene\t{g.start}\t{g.end}\t.\t{g.strand}\t.\t"
                 f"{_fmt_gff3_attrs({'ID': g.id, **g.attrs})}\n")
        for t in sorted(g.transcripts.values(), key=lambda t: (t.start, t.id)):
            ttype = "mRNA" if t.is_coding else "transcript"
            fh.write(f"{t.contig}\t{source}\t{ttype}\t{t.start}\t{t.end}\t.\t{t.strand}\t.\t"
                     f"{_fmt_gff3_attrs({'ID': t.id, 'Parent': g.id, **t.attrs})}\n")
            five, three = t.utrs()
            for i, (a, b) in enumerate(t.exons, 1):
                fh.write(f"{t.contig}\t{source}\texon\t{a}\t{b}\t.\t{t.strand}\t.\t"
                         f"ID={quote(t.id)}.exon{i};Parent={quote(t.id)}\n")
            for kind, spans in (("five_prime_UTR", five), ("three_prime_UTR", three)):
                for a, b in spans:
                    fh.write(f"{t.contig}\t{source}\t{kind}\t{a}\t{b}\t.\t{t.strand}\t.\t"
                             f"Parent={quote(t.id)}\n")
            for (a, b), ph in zip(t.cds, t.cds_phases()):
                fh.write(f"{t.contig}\t{source}\tCDS\t{a}\t{b}\t.\t{t.strand}\t{ph}\t"
                         f"ID=cds.{quote(t.id)};Parent={quote(t.id)}\n")
        fh.write("\n")


def write_gtf(genes: Iterable[Gene], fh: TextIO, source: str = "PASA-XG") -> None:
    for g in sorted(genes, key=lambda g: (g.contig, g.start, g.id)):
        for t in sorted(g.transcripts.values(), key=lambda t: (t.start, t.id)):
            ids = f'gene_id "{g.id}"; transcript_id "{t.id}";'
            fh.write(f"{t.contig}\t{source}\ttranscript\t{t.start}\t{t.end}\t.\t{t.strand}\t.\t{ids}\n")
            for a, b in t.exons:
                fh.write(f"{t.contig}\t{source}\texon\t{a}\t{b}\t.\t{t.strand}\t.\t{ids}\n")
            for (a, b), ph in zip(t.cds, t.cds_phases()):
                fh.write(f"{t.contig}\t{source}\tCDS\t{a}\t{b}\t.\t{t.strand}\t{ph}\t{ids}\n")


def write_bed(transcripts: Iterable[Transcript], fh: TextIO) -> None:
    for t in sorted(transcripts, key=lambda t: (t.contig, t.start, t.id)):
        thick_s, thick_e = (t.cds[0][0] - 1, t.cds[-1][1]) if t.cds else (t.start - 1, t.start - 1)
        sizes = ",".join(str(b - a + 1) for a, b in t.exons)
        starts = ",".join(str(a - t.start) for a, b in t.exons)
        fh.write(f"{t.contig}\t{t.start - 1}\t{t.end}\t{t.id}\t0\t{t.strand}\t{thick_s}\t{thick_e}\t"
                 f"0\t{t.num_exons}\t{sizes}\t{starts}\n")


def write_fasta(records: Iterable[Tuple[str, str]], fh: TextIO, width: int = 60) -> None:
    for name, seq in records:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")


def read_fasta(path: str) -> Dict[str, str]:
    out: Dict[str, List[str]] = OrderedDict()
    name: Optional[str] = None
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = []
            elif name is not None:
                out[name].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}
