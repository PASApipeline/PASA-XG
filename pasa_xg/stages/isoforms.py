"""Stage 3: load LRAA's isoforms into one table, write their sequences."""

from __future__ import annotations

import os
from typing import Any, Dict

from .. import isoforms as iso_mod
from ..config import Settings
from ..genome import Genome
from ..gff import write_fasta
from ..runner import Run, code_digest


def run(run: Run, s: Settings, prep: Dict[str, Any], asm: Dict[str, Any]) -> Dict[str, Any]:
    inputs = {"gtf": asm["gtf_path"], "expr": asm["expr_path"], "genome": prep["genome_path"],
              "stage": "isoforms/1", "src": code_digest("pasa_xg.stages.isoforms", "pasa_xg.isoforms", "pasa_xg.gff")}

    def body(d: str) -> Dict[str, Any]:
        isos = iso_mod.load_lraa(asm["gtf_path"], asm["expr_path"])
        genome = Genome(prep["genome_path"])
        fa = os.path.join(d, "isoforms.fa")
        with open(fa, "w") as fh:
            write_fasta(((i.id, genome.transcript_seq(i.tx)) for i in isos.values()), fh)
        pkl = os.path.join(d, "isoforms.pkl")
        iso_mod.save(isos, pkl)
        iso_mod.table(isos).to_parquet(os.path.join(d, "isoforms.parquet"), index=False)
        return {"isoforms_path": pkl, "fasta_path": fa,
                "parquet_path": os.path.join(d, "isoforms.parquet"),
                "isoforms": len(isos),
                "multi_exon": sum(1 for i in isos.values() if i.tx.num_exons > 1)}

    return run.run_stage("isoforms", "03_isoforms", inputs, body)
