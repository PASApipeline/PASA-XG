"""Stage 8: collect the final deliverables in outputs/ and write a summary table."""

from __future__ import annotations

import os
import shutil
from typing import Any, Dict

import pandas as pd

from ..runner import Run

DELIVERABLES = [
    ("03_isoforms", "isoforms.fa", "isoforms.fa"),
    ("04_coding", "isoforms.coding.parquet", "isoforms.parquet"),
    ("05_annot_update", "annotation.updated.gff3", "annotation.updated.gff3"),
    ("05_annot_update", "annotation.updated.gtf", "annotation.updated.gtf"),
    ("05_annot_update", "annotation.updated.bed", "annotation.updated.bed"),
    ("05_annot_update", "annotation.updated.pep", "annotation.updated.pep"),
    ("05_annot_update", "annot_update.decisions.tsv", "annot_update.decisions.tsv"),
    ("05_annot_update", "annot_update.rounds.tsv", "annot_update.rounds.tsv"),
    ("06_gene_models", "gene_models.gff3", "gene_models.gff3"),
    ("06_gene_models", "gene_models.pep", "gene_models.pep"),
]


def run(run: Run, results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    out_dir = run.path("outputs")
    os.makedirs(out_dir, exist_ok=True)
    for sub, name, dest in DELIVERABLES:
        src = run.path(sub, name)
        if os.path.exists(src):
            shutil.copyfile(src, os.path.join(out_dir, dest))
    asm = results.get("assemble", {})
    for k in ("gtf_path", "expr_path", "tracking_path"):
        if k in asm and os.path.exists(asm[k]):
            shutil.copyfile(asm[k], os.path.join(out_dir, os.path.basename(asm[k])))
    rows = [{"stage": stage, "metric": k, "value": v}
            for stage, o in results.items() for k, v in o.items()
            if not k.endswith("_path") and isinstance(v, (int, float, str, bool))]
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, "summary.tsv"), sep="\t", index=False)
    run.manifest["stages"]["report"] = {"status": "done"}
    run.save_manifest()
    return {"outputs_path": out_dir}
