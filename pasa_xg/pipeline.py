"""Stage order and the full run."""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional

from .config import Settings
from .runner import Run
from .stages import align, assemble, coding, gene_models, isoforms, prep

log = logging.getLogger("pasa_xg")

STAGES = ["prep", "align", "assemble", "isoforms", "coding", "annot_update", "gene_models"]


def run_all(s: Settings, force_from: Optional[str] = None, stop_after: Optional[str] = None) -> Dict[str, Any]:
    if force_from and force_from not in STAGES:
        raise ValueError(f"unknown stage {force_from}; stages: {', '.join(STAGES)}")
    run = Run(s.out_dir, s.to_dict(), force_from=force_from)
    out: Dict[str, Any] = {}

    def done(name: str) -> bool:
        return stop_after == name

    out["prep"] = prep.run(run, s)
    if done("prep"):
        return out
    out["align"] = align.run(run, s, out["prep"])
    if done("align"):
        return out
    out["assemble"] = assemble.run(run, s, out["prep"], out["align"])
    if done("assemble"):
        return out
    out["isoforms"] = isoforms.run(run, s, out["prep"], out["assemble"])
    if done("isoforms"):
        return out
    out["coding"] = coding.run(run, s, out["isoforms"])
    if done("coding"):
        return out
    if s.annotation and out["prep"].get("coding_transcripts", 0) > 0:
        from .stages import annot_update
        out["annot_update"] = annot_update.run(run, s, out["prep"], out["coding"])
    elif s.annotation:
        log.warning("annotation has no CDS features: annotation update skipped")
    if done("annot_update"):
        return out
    out["gene_models"] = gene_models.run(run, s, out["prep"], out["coding"])
    return out
