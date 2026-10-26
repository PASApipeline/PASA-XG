"""Stage 1: align reads or transcripts with minimap2 (skipped for BAM input)."""

from __future__ import annotations

import os
import shlex
from typing import Any, Dict

from ..config import Settings
from ..runner import Run, run_cmd


def run(run: Run, s: Settings, prep: Dict[str, Any]) -> Dict[str, Any]:
    inputs = {"mmi": prep.get("mmi_path"), "genome": prep["genome_path"],
              "reads": [run.fdigest(r) for r in s.reads], "bam": s.bam,
              "args": s.minimap2_args(), "G": s.max_intron_length, "stage": "align/1"}

    def body(d: str) -> Dict[str, Any]:
        bam = os.path.join(d, "aligned.bam")
        logp = os.path.join(d, "align.log")
        if s.bam:
            if os.path.lexists(bam):
                os.remove(bam)
            os.symlink(os.path.abspath(s.bam), bam)
            if os.path.exists(os.path.abspath(s.bam) + ".bai"):
                os.symlink(os.path.abspath(s.bam) + ".bai", bam + ".bai")
            else:
                run_cmd(["samtools", "index", bam], log_path=logp)
            return {"bam_path": bam, "aligned": False}
        mm2 = ["minimap2", *s.minimap2_args(), "-G", str(s.max_intron_length), "-t", str(s.cpu),
               prep["mmi_path"], *s.reads]
        sort = ["samtools", "sort", "-@", str(max(1, s.cpu // 2)), "-o", bam, "-"]
        cmd = ("set -o pipefail; " + " ".join(shlex.quote(c) for c in mm2)
               + f" 2>> {shlex.quote(logp)} | " + " ".join(shlex.quote(c) for c in sort))
        run_cmd(cmd, log_path=logp)
        run_cmd(["samtools", "index", bam], log_path=logp)
        return {"bam_path": bam, "aligned": True}

    return run.run_stage("align", "01_align", inputs, body)
