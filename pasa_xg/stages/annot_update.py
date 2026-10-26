"""Stage 5: compare the isoforms with the annotation and update it, in rounds.

LRAA's isoforms are fixed. Round 1 compares every isoform; each later round
revisits only the LRAA genes that overlap a gene the previous round changed or
created and that still hold unincorporated isoforms. Stops after --max_rounds or
when a round incorporates nothing new.
"""

from __future__ import annotations

import logging
import os
import pickle
from collections import Counter
from typing import Any, Dict, List, Set

import pandas as pd

from .. import isoforms as iso_mod
from ..annot.apply import apply_round
from ..annot.checks import Params
from ..annot.compare import STATUS, Decision, run_round
from ..config import Settings
from ..genome import Genome, codon_table, translate
from ..gff import write_bed, write_fasta, write_gff3, write_gtf
from ..models import AnnotationSet
from ..runner import Run, code_digest

log = logging.getLogger("pasa_xg")


def _rows(decisions: List[Decision], applied: Dict[str, str], rnd: int, isos) -> List[Dict[str, Any]]:
    out = []
    for d in decisions:
        i = isos[d.isoform_id]
        out.append({"round": rnd, "isoform_id": d.isoform_id, "lraa_gene": d.group,
                    "action": d.action, "valid": d.valid, "pasa_status": d.status,
                    "pasa_status_descr": STATUS.get(d.status, ""), "genes": " ".join(d.genes),
                    "model": d.model or "", "annotation_id": applied.get(d.isoform_id, ""),
                    "full_length": i.full_length, "num_exons": i.tx.num_exons,
                    "TPM": i.tpm, "isoform_fraction": i.isoform_fraction, "comment": d.comment})
    return out


def write_annotation(ann: AnnotationSet, genome: Genome, table, prefix: str) -> None:
    with open(prefix + ".gff3", "w") as fh:
        write_gff3(ann, fh)
    with open(prefix + ".gtf", "w") as fh:
        write_gtf(ann, fh)
    with open(prefix + ".bed", "w") as fh:
        write_bed(ann.transcripts(), fh)
    coding = [t for t in ann.transcripts() if t.cds]
    with open(prefix + ".pep", "w") as fh:
        write_fasta(((t.id, genome.protein(t, table)) for t in coding), fh)


def update(ann: AnnotationSet, isos: Dict[str, iso_mod.Isoform], genome: Genome, s: Settings,
           out_dir: str, p: Params) -> Dict[str, Any]:
    table = codon_table(s.genetic_code)
    expression = {k: {"TPM": v.tpm} for k, v in isos.items()}
    eligible = {k: v for k, v in isos.items()
                if v.tpm >= s.min_tpm and v.isoform_fraction >= s.min_isoform_fraction}
    rows: List[Dict[str, Any]] = []
    for k in set(isos) - set(eligible):
        i = isos[k]
        rows.append({"round": 0, "isoform_id": k, "lraa_gene": i.tx.gene_id, "action": "below_threshold",
                     "valid": False, "pasa_status": 100, "pasa_status_descr": STATUS[100],
                     "genes": "", "model": "", "annotation_id": "", "full_length": i.full_length,
                     "num_exons": i.tx.num_exons, "TPM": i.tpm,
                     "isoform_fraction": i.isoform_fraction, "comment": ""})
    incorporated: Set[str] = set()
    support: Dict[str, Set[str]] = {}   # annotation transcript -> isoforms that built/confirm it
    fl_models: Set[str] = set()         # ... of which a full-length isoform
    targets = set(eligible)
    summary = []
    current = ann
    for rnd in range(1, s.max_rounds + 1):
        if not targets:
            break
        rdir = os.path.join(out_dir, "rounds", f"round_{rnd:02d}")
        os.makedirs(rdir, exist_ok=True)
        decisions, splits, _ = run_round(current, [eligible[k] for k in sorted(targets)],
                                         genome, table, p, support, fl_models,
                                         all_isoforms=eligible.values())
        label = f"round_{rnd:02d}"
        new, applied, changed = apply_round(current, decisions, splits, label, expression)
        newly = {d.isoform_id for d in decisions
                 if d.valid and d.action not in ("deferred", "antisense", "intergenic_noncoding")
                 and d.isoform_id in applied}
        incorporated |= newly
        for iso_id in newly:
            tid = applied[iso_id]
            support.setdefault(tid, set()).add(iso_id)
            if isos[iso_id].full_length:
                fl_models.add(tid)
        rrows = _rows(decisions, applied, rnd, isos)
        rows += rrows
        pd.DataFrame(rrows).to_csv(os.path.join(rdir, "decisions.tsv"), sep="\t", index=False)
        write_annotation(new, genome, table, os.path.join(rdir, "annotation"))
        counts = Counter(r["action"] for r in rrows)
        summary.append({"round": rnd, "evaluated": len(decisions), "incorporated": len(newly),
                        "genes_changed": len(changed), **counts})
        log.info("annotation update round %d: %d isoforms evaluated, %d incorporated, %d genes changed",
                 rnd, len(decisions), len(newly), len(changed))
        current = new
        if not newly:
            break
        # next round: unincorporated isoforms of LRAA genes overlapping changed genes
        groups = {}
        for k in eligible:
            if k not in incorporated:
                groups.setdefault(eligible[k].tx.gene_id, []).append(k)
        targets = set()
        for gid in changed:
            if gid not in current.genes:
                continue
            g = current.genes[gid]
            for grp, ks in groups.items():
                t0 = eligible[ks[0]].tx
                if t0.contig != g.contig:
                    continue
                lo = min(eligible[k].tx.start for k in ks)
                hi = max(eligible[k].tx.end for k in ks)
                if lo <= g.end and hi >= g.start:
                    targets |= set(ks)
    # final state per isoform: its last evaluation
    table_df = pd.DataFrame(rows)
    final = table_df.sort_values("round").groupby("isoform_id", as_index=False).last()
    final.loc[final["isoform_id"].isin(incorporated), "incorporated"] = True
    final["incorporated"] = final["isoform_id"].isin(incorporated)
    final.to_csv(os.path.join(out_dir, "annot_update.decisions.tsv"), sep="\t", index=False)
    table_df.to_csv(os.path.join(out_dir, "annot_update.all_rounds.tsv"), sep="\t", index=False)
    pd.DataFrame(summary).to_csv(os.path.join(out_dir, "annot_update.rounds.tsv"), sep="\t", index=False)
    write_annotation(current, genome, table, os.path.join(out_dir, "annotation.updated"))
    with open(os.path.join(out_dir, "annotation.updated.pkl"), "wb") as fh:
        pickle.dump(current, fh)
    return {"rounds": len(summary), "incorporated": len(incorporated),
            "gff3_path": os.path.join(out_dir, "annotation.updated.gff3")}


def run(run: Run, s: Settings, prep: Dict[str, Any], coding: Dict[str, Any]) -> Dict[str, Any]:
    inputs = {"annotation": prep["annotation_models_path"], "isoforms": coding["isoforms_path"],
              "rounds": s.max_rounds, "protect": s.protect_supported_models, "min_tpm": s.min_tpm, "min_frac": s.min_isoform_fraction,
              "code": s.genetic_code, "stage": "annot_update/1", "src": code_digest("pasa_xg.stages.annot_update", "pasa_xg.annot", "pasa_xg.models", "pasa_xg.genome", "pasa_xg.gff")}

    def body(d: str) -> Dict[str, Any]:
        with open(prep["annotation_models_path"], "rb") as fh:
            ann = pickle.load(fh)
        isos = iso_mod.load(coding["isoforms_path"])
        return update(ann, isos, Genome(prep["genome_path"]), s, d,
                      Params(PROTECT_SUPPORTED_MODELS=s.protect_supported_models))

    return run.run_stage("annot_update", "05_annot_update", inputs, body)
