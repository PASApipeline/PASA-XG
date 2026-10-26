"""Stage 4: TransDecoder ORFs on the isoforms, mapped to the genome."""

from __future__ import annotations

import os
import re
from typing import Any, Dict, Optional

from .. import isoforms as iso_mod
from ..config import Settings
from ..isoforms import Orf
from ..gff import read_fasta
from ..runner import Run, code_digest, run_cmd

# >tid.p1 GENE.x~~tid.p1  ORF type:complete (+),score=404.17 len:483 tid:410-1861(+)
_PEP_HEADER = re.compile(
    r"^>(\S+)\s.*ORF type:(\S+) \(([+-])\),score=(\S+)\s+len:(\d+)\s+(\S+):(\d+)-(\d+)\(([+-])\)\s*$")

# TransDecoder's start-codon refinement trains a model on the longest ORFs and
# fails outright on small inputs, so it is skipped below this many isoforms.
MIN_ISOFORMS_FOR_START_REFINEMENT = 1000


def parse_pep_headers(pep_path: str) -> Dict[str, Orf]:
    """Best ORF per transcript from TransDecoder's .pep headers."""
    best: Dict[str, Orf] = {}
    with open(pep_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            m = _PEP_HEADER.match(line.rstrip("\n"))
            if not m:
                raise ValueError(f"unrecognized TransDecoder header: {line.strip()}")
            orf_id, otype, ostrand, score, aa_len, tid, a, b, _ = m.groups()
            if ostrand != "+":
                continue  # -S: isoforms are oriented, minus-strand ORFs are not wanted
            orf = Orf(orf_id, otype, int(a), int(b), int(aa_len), float(score))
            if tid not in best or orf.aa_len > best[tid].aa_len:
                best[tid] = orf
    return best


def refine_start(orf: Orf, seq: str, min_aa: int) -> Orf:
    """A 5'-partial (or internal) ORF on a full-length isoform usually means TransDecoder
    ran the ORF off the 5' end past its real start codon (its own start refinement
    needs a large training set). Start it at the first in-frame ATG instead, when the
    protein keeps at least `min_aa` residues."""
    if orf.orf_type not in ("5prime_partial", "internal"):
        return orf
    for pos in range(orf.t_start, orf.t_end - 2, 3):  # 1-based codon starts
        if seq[pos - 1:pos + 2].upper() == "ATG":
            aa = (orf.t_end - pos + 1) // 3
            if aa < min_aa:
                return orf
            kind = "complete" if orf.orf_type == "5prime_partial" else "3prime_partial"
            return Orf(orf.orf_id, kind, pos, orf.t_end, aa - (1 if kind == "complete" else 0),
                       orf.score)
    return orf


def transdecoder_bin(s: Settings, name: str) -> str:
    return os.path.join(s.transdecoder_dir, name) if s.transdecoder_dir else name


def run(run: Run, s: Settings, isos_out: Dict[str, Any]) -> Dict[str, Any]:
    inputs = {"isoforms": isos_out["isoforms_path"], "min_aa": s.min_orf_aa,
              "code": s.genetic_code, "pfam": s.pfam_hits, "blastp": s.blastp_hits,
              "stage": "coding/1", "src": code_digest("pasa_xg.stages.coding", "pasa_xg.isoforms", "pasa_xg.models")}

    def body(d: str) -> Dict[str, Any]:
        logp = os.path.join(d, "transdecoder.log")
        fa = os.path.join(d, "isoforms.fa")
        if os.path.lexists(fa):
            os.remove(fa)
        os.symlink(isos_out["fasta_path"], fa)
        isos = iso_mod.load(isos_out["isoforms_path"])
        td = os.path.join(d, "transdecoder")
        run_cmd([transdecoder_bin(s, "TransDecoder.LongOrfs"), "-t", fa, "-S", "-O", td,
                 "-m", str(s.min_orf_aa), "-G", s.genetic_code], log_path=logp, cwd=d)
        predict = [transdecoder_bin(s, "TransDecoder.Predict"), "-t", fa, "-O", td,
                   "--single_best_only", "-G", s.genetic_code]
        if len(isos) < MIN_ISOFORMS_FOR_START_REFINEMENT:
            predict.append("--no_refine_starts")
        if s.pfam_hits:
            predict += ["--retain_pfam_hits", s.pfam_hits]
        if s.blastp_hits:
            predict += ["--retain_blastp_hits", s.blastp_hits]
        run_cmd(predict, log_path=logp, cwd=d)
        orfs = parse_pep_headers(os.path.join(td, "isoforms.fa.transdecoder.pep"))
        seqs = read_fasta(isos_out["fasta_path"])
        refined = 0
        for tid in list(orfs):
            new = refine_start(orfs[tid], seqs[tid], s.min_orf_aa)
            refined += new is not orfs[tid]
            orfs[tid] = new
        for tid, iso in isos.items():
            orf: Optional[Orf] = orfs.get(tid)
            iso.orf = orf
            iso.tx.cds = iso.tx.segments_for(orf.t_start, orf.t_end) if orf else []
        pkl = os.path.join(d, "isoforms.coding.pkl")
        iso_mod.save(isos, pkl)
        iso_mod.table(isos).to_parquet(os.path.join(d, "isoforms.coding.parquet"), index=False)
        n_complete = sum(1 for i in isos.values() if i.full_length)
        return {"isoforms_path": pkl, "with_orf": len(orfs), "complete_orfs": n_complete,
                "starts_refined": refined}

    return run.run_stage("coding", "04_coding", inputs, body)
