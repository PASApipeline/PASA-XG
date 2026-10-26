"""Stage 0: index the genome, load and check the annotation."""

from __future__ import annotations

import logging
import os
import pickle
from typing import Any, Dict

import pandas as pd

from ..annot_check import check_annotation
from ..config import Settings
from ..genome import Genome
from ..gff import read_annotation, write_gtf
from ..runner import Run, run_cmd

log = logging.getLogger("pasa_xg")


def link_genome(s: Settings, d: str) -> str:
    """Symlink the genome into the stage dir and faidx it there."""
    fa = os.path.join(d, "genome.fa")
    if os.path.lexists(fa):
        os.remove(fa)
    os.symlink(os.path.abspath(s.genome), fa)
    for ext in (".fai",):
        if os.path.exists(fa + ext):
            os.remove(fa + ext)
    run_cmd(["samtools", "faidx", fa], log_path=os.path.join(d, "prep.log"))
    return fa


def run(run: Run, s: Settings) -> Dict[str, Any]:
    inputs = {"genome": s.genome, "annotation": s.annotation, "code": s.genetic_code,
              "align": bool(s.reads), "stage": "prep/1"}

    def body(d: str) -> Dict[str, Any]:
        out: Dict[str, Any] = {"genome_path": link_genome(s, d)}
        if s.reads:
            mmi = os.path.join(d, "genome.mmi")
            run_cmd(["minimap2", "-d", mmi, "-t", str(s.cpu), out["genome_path"]],
                    log_path=os.path.join(d, "prep.log"))
            out["mmi_path"] = mmi
        if s.annotation:
            ann = read_annotation(s.annotation)
            genome = Genome(out["genome_path"])
            missing = {g.contig for g in ann} - set(genome.contigs)
            if missing:
                raise ValueError(f"annotation contigs absent from the genome: {sorted(missing)[:5]}")
            checks = check_annotation(ann, genome, s.genetic_code)
            table = pd.DataFrame([c.row() for c in checks])
            table.to_csv(os.path.join(d, "annotation_check.tsv"), sep="\t", index=False)
            n_coding = int(table["coding"].sum()) if len(table) else 0
            n_bad = int((~table["ok"]).sum()) if len(table) else 0
            with open(os.path.join(d, "annotation.models.pkl"), "wb") as fh:
                pickle.dump(ann, fh)
            with open(os.path.join(d, "annotation.gtf"), "w") as fh:
                write_gtf(ann, fh)
            log.info("annotation: %d genes, %d transcripts, %d coding, %d with CDS problems",
                     len(ann), len(table), n_coding, n_bad)
            out.update(annotation_models_path=os.path.join(d, "annotation.models.pkl"),
                       annotation_gtf_path=os.path.join(d, "annotation.gtf"),
                       annotation_check_path=os.path.join(d, "annotation_check.tsv"),
                       genes=len(ann), transcripts=len(table), coding_transcripts=n_coding,
                       cds_problems=n_bad)
        return out

    return run.run_stage("prep", "00_prep", inputs, body)
