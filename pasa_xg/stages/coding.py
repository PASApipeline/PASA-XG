"""Stage 4: TransDecoder ORFs on the isoforms (pytransdecoder, the Python TransDecoder), mapped to the genome."""

from __future__ import annotations

import os
import re
from typing import Any, Dict, Optional

from .. import isoforms as iso_mod
from ..config import Settings
from ..isoforms import Orf
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


def diamond_hits(s: Settings, pep: str, d: str, logp: str) -> str:
    """diamond blastp of the long ORFs against --homology_proteins (outfmt 6, best hit per ORF, E <= 1e-5)."""
    db = os.path.join(d, "homology_proteins")
    out = os.path.join(d, "longest_orfs.diamond.outfmt6")
    run_cmd([s.diamond_cmd, "makedb", "--in", os.path.abspath(s.homology_proteins), "-d", db, "--threads", str(s.cpu)],
            log_path=logp, cwd=d)
    run_cmd([s.diamond_cmd, "blastp", "-q", pep, "-d", db, "-o", out, "--outfmt", "6", "--evalue", "1e-5",
             "--max-target-seqs", "1", "--threads", str(s.cpu)], log_path=logp, cwd=d)
    return out


def run(run: Run, s: Settings, isos_out: Dict[str, Any]) -> Dict[str, Any]:
    inputs = {"isoforms": isos_out["isoforms_path"], "min_aa": s.min_orf_aa,
              "code": s.genetic_code, "pfam": s.pfam_hits, "blastp": s.blastp_hits,
              "homology": run.fdigest(s.homology_proteins) if s.homology_proteins else None,
              "stage": "coding/4", "transdecoder": s.transdecoder_cmd, "src": code_digest("pasa_xg.stages.coding", "pasa_xg.isoforms", "pasa_xg.models")}

    def body(d: str) -> Dict[str, Any]:
        logp = os.path.join(d, "transdecoder.log")
        fa = os.path.join(d, "isoforms.fa")
        if os.path.lexists(fa):
            os.remove(fa)
        os.symlink(isos_out["fasta_path"], fa)
        isos = iso_mod.load(isos_out["isoforms_path"])
        td = os.path.join(d, "transdecoder")
        run_cmd([s.transdecoder_cmd, "longorfs", "-t", fa, "-S", "-O", td,
                 "-m", str(s.min_orf_aa), "-G", s.genetic_code], log_path=logp, cwd=d)
        # one ORF per transcript is pytransdecoder's default (Perl's --single_best_only)
        predict = [s.transdecoder_cmd, "predict", "-t", fa, "-O", td, "-G", s.genetic_code]
        if len(isos) < MIN_ISOFORMS_FOR_START_REFINEMENT:
            predict.append("--no_refine_starts")
        if s.pfam_hits:
            predict += ["--retain_pfam_hits", s.pfam_hits]
        blastp = [s.blastp_hits] if s.blastp_hits else []
        if s.homology_proteins:
            blastp.append(diamond_hits(s, os.path.join(td, "isoforms.fa.transdecoder_dir", "longest_orfs.pep"),
                                       d, logp))
        if len(blastp) > 1:
            both = os.path.join(d, "blastp_hits.outfmt6")
            with open(both, "w") as out:
                for path in blastp:
                    with open(path) as fh:
                        out.write(fh.read())
            blastp = [both]
        if blastp:
            predict += ["--retain_blastp_hits", blastp[0]]
        run_cmd(predict, log_path=logp, cwd=d)
        orfs = parse_pep_headers(os.path.join(td, "isoforms.fa.transdecoder.pep"))
        for tid, iso in isos.items():
            orf: Optional[Orf] = orfs.get(tid)
            iso.orf = orf
            iso.tx.cds = iso.tx.segments_for(orf.t_start, orf.t_end) if orf else []
        pkl = os.path.join(d, "isoforms.coding.pkl")
        iso_mod.save(isos, pkl)
        iso_mod.table(isos).to_parquet(os.path.join(d, "isoforms.coding.parquet"), index=False)
        n_complete = sum(1 for i in isos.values() if i.full_length)
        return {"isoforms_path": pkl, "with_orf": len(orfs), "complete_orfs": n_complete}

    return run.run_stage("coding", "04_coding", inputs, body)
