"""Stage 2: assemble and quantify isoforms with LRAA."""

from __future__ import annotations

import glob
import os
import shlex
from typing import Any, Dict

from ..config import Settings
from ..runner import Run, run_cmd, tool_version


def _one(d: str, pattern: str) -> str:
    hits = [h for h in glob.glob(os.path.join(d, pattern)) if "splice_pattern_collapsed" not in h]
    if len(hits) != 1:
        raise RuntimeError(f"expected one LRAA output matching {pattern} in {d}, found {hits}")
    return hits[0]


def run(run: Run, s: Settings, prep: Dict[str, Any], align: Dict[str, Any]) -> Dict[str, Any]:
    lraa = shlex.split(s.lraa_cmd)
    version = tool_version([*lraa, "--version"])
    inputs = {"bam": align["bam_path"], "genome": prep["genome_path"],
              "gtf": prep.get("annotation_gtf_path"), "args": s.lraa_args(),
              "lraa": version, "stage": "assemble/1"}

    def body(d: str) -> Dict[str, Any]:
        cmd = [*lraa, "--genome", prep["genome_path"], "--bam", align["bam_path"],
               "--output_prefix", "lraa", "--cpu_budget", str(s.cpu), *s.lraa_args()]
        if prep.get("annotation_gtf_path"):
            cmd += ["--gtf", prep["annotation_gtf_path"]]
        run_cmd(cmd, log_path=os.path.join(d, "lraa.log"), cwd=d)
        return {"gtf_path": _one(d, "lraa.LRAA.*.gtf"),
                "expr_path": _one(d, "lraa.LRAA.*.quant.expr"),
                "tracking_path": _one(d, "lraa.LRAA.*.quant.tracking*"),
                "lraa_version": version}

    return run.run_stage("assemble", "02_assemble", inputs, body)
