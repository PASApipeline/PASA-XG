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
from ..annot import refframe
from ..annot.noncoding import NC_STATUS, run_noncoding, tie_kind
from ..annot.structure import introns_of
from ..annot.sites import POLYA, TSS, SiteIndex, add_isoform_ends, five_prime, load_bed, three_prime
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


def site_index(s: Settings, isoforms) -> SiteIndex:
    """Known sites from --tss_sites / --polya_sites plus the isoforms' LRAA-called ends."""
    idx = SiteIndex(s.site_window)
    for paths, kind in ((s.tss_sites, TSS), (s.polya_sites, POLYA)):
        for path in paths:
            log.info("%s sites: %d from %s", kind, load_bed(idx, path, kind), path)
    log.info("LRAA-called isoform ends: %d", add_isoform_ends(idx, isoforms))
    return idx


def label_ends(ann: AnnotationSet, idx: SiteIndex) -> None:
    """tss_support / polya_support on every transcript: the source of the closest site
    within the window of its 5' / 3' end ('lraa', 'known:<bed>'), or 'none'."""
    for t in ann.transcripts():
        t.attrs["tss_support"] = idx.near(t.contig, t.strand, TSS, five_prime(t)) or "none"
        t.attrs["polya_support"] = idx.near(t.contig, t.strand, POLYA, three_prime(t)) or "none"


TIE_STATUS = {5, 11, 30, 35, 105}   # matches / overlaps two or more genes; not updated


def write_conflicts(rows: List[Dict[str, Any]], incorporated: Set[str], isos, ann: AnnotationSet,
                    path: str) -> None:
    """Gene sets that unincorporated isoforms tie between, in the coding rounds (5, 11, 30, 35) or the
    noncoding pass (105), for the user to review: the annotation may split one gene
    (shared_introns, or genes_share_introns), or the isoforms read through (readthrough) or only
    overlap the genes' exons (overlap)."""
    last: Dict[str, Dict[str, Any]] = {}
    for r in rows:
        if r["pasa_status"] in TIE_STATUS and r["isoform_id"] not in incorporated:
            last[r["isoform_id"]] = r
    by: Dict[Any, List[Dict[str, Any]]] = {}
    for k, r in last.items():
        genes = [ann.genes[g] for g in sorted(set(r["genes"].split())) if g in ann.genes]
        if len(genes) < 2:
            continue
        by.setdefault((" ".join(g.id for g in genes), tie_kind(isos[k].tx, genes)), []).append(r)
    def share(gs: str) -> bool:
        per = [set().union(*(introns_of(t.exons) for t in ann.genes[g].transcripts.values())) for g in gs.split()]
        return bool(set.intersection(*per))

    out = [{"genes": g, "kind": kind, "genes_share_introns": share(g), "isoforms": len(rs),
            "pasa_status": ",".join(str(x) for x in sorted({r["pasa_status"] for r in rs})),
            "max_TPM": max(isos[r["isoform_id"]].tpm for r in rs),
            "isoform_ids": ",".join(r["isoform_id"] for r in rs)} for (g, kind), rs in by.items()]
    df = pd.DataFrame(out, columns=["genes", "kind", "genes_share_introns", "isoforms", "pasa_status", "max_TPM", "isoform_ids"])
    df.sort_values(["kind", "isoforms"], ascending=[False, False]).to_csv(path, sep="\t", index=False)
    if len(df):
        log.warning("%d unincorporated isoforms tie between two or more genes (%s gene sets); not updated, "
                    "see %s", int(df["isoforms"].sum()),
                    ", ".join(f"{n} {k}" for k, n in df["kind"].value_counts().items()), os.path.basename(path))


def update(ann: AnnotationSet, isos: Dict[str, iso_mod.Isoform], genome: Genome, s: Settings,
           out_dir: str, p: Params) -> Dict[str, Any]:
    table = codon_table(s.genetic_code)
    expression = {k: {"TPM": v.tpm} for k, v in isos.items()}
    eligible = {k: v for k, v in isos.items()
                if v.tpm >= s.min_tpm and v.isoform_fraction >= s.min_isoform_fraction}
    sites = site_index(s, eligible.values())
    if s.reference_frame_check:
        refframe.check(eligible, ann, genome, table, sites if len(sites) else None, p, s.min_orf_aa)
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
        decisions, splits, rc = run_round(current, [eligible[k] for k in sorted(targets)],
                                          genome, table, p, support, fl_models,
                                          all_isoforms=eligible.values(), sites=sites)
        label = f"round_{rnd:02d}"
        new, applied, changed = apply_round(current, decisions, splits, label, expression, s.run_label,
                                            alt_sites=rc.alt_sites)
        newly = {d.isoform_id for d in decisions
                 if d.valid and d.action not in ("deferred", "antisense", "intergenic_noncoding")
                 and d.action != "alt_site_recorded" and d.isoform_id in applied}
        recorded = {d.isoform_id for d in decisions if d.valid and d.action == "alt_site_recorded"}
        incorporated |= recorded   # evidence kept on the model; not revisited
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
    if s.include_noncoding:
        last = {}
        for r in rows:
            last[r["isoform_id"]] = r
        cands = [k for k in eligible if k not in incorporated and last.get(k, {}).get("action") != "deferred"]
        current, ncd = run_noncoding(current, eligible, cands, genome, table, p, sites if len(sites) else None,
                                     s.run_label)
        rnd = len(summary) + 1
        nrows = []
        for d in ncd:
            i, prev = isos[d.isoform_id], last.get(d.isoform_id, {})
            nrows.append({"round": rnd, "isoform_id": d.isoform_id, "lraa_gene": i.tx.gene_id,
                          "action": d.action, "valid": d.status in (101, 102), "pasa_status": d.status,
                          "pasa_status_descr": NC_STATUS[d.status], "genes": " ".join(d.genes), "model": "",
                          "annotation_id": d.annotation_id, "full_length": i.full_length,
                          "num_exons": i.tx.num_exons, "TPM": i.tpm, "isoform_fraction": i.isoform_fraction,
                          "comment": "; ".join(x for x in (d.comment, f"coding round: {prev.get('action', '')} "
                                                           f"({prev.get('pasa_status', '')}) {prev.get('comment', '')}".strip()) if x)})
            if d.status in (101, 102):
                incorporated.add(d.isoform_id)
        rows += nrows
        counts = Counter(r["action"] for r in nrows)
        summary.append({"round": rnd, "evaluated": len(nrows), "incorporated": counts["noncoding_isoform"]
                        + counts["noncoding_gene"], "genes_changed": len({r["genes"] for r in nrows if r["valid"]}),
                        **counts})
    if s.reference_frame_check:
        refframe.check_models(current, ann, genome, table, sites if len(sites) else None, p, s.min_orf_aa)
    refframe.sync_partial_tags(current, ann, genome, table)
    write_conflicts(rows, incorporated, isos, current, os.path.join(out_dir, "annot_update.gene_conflicts.tsv"))
    # final state per isoform: its last evaluation
    table_df = pd.DataFrame(rows)
    final = table_df.sort_values("round").groupby("isoform_id", as_index=False).last()
    final.loc[final["isoform_id"].isin(incorporated), "incorporated"] = True
    final["incorporated"] = final["isoform_id"].isin(incorporated)
    final.to_csv(os.path.join(out_dir, "annot_update.decisions.tsv"), sep="\t", index=False)
    table_df.to_csv(os.path.join(out_dir, "annot_update.all_rounds.tsv"), sep="\t", index=False)
    pd.DataFrame(summary).to_csv(os.path.join(out_dir, "annot_update.rounds.tsv"), sep="\t", index=False)
    if len(sites):
        label_ends(current, sites)
    write_annotation(current, genome, table, os.path.join(out_dir, "annotation.updated"))
    with open(os.path.join(out_dir, "annotation.updated.pkl"), "wb") as fh:
        pickle.dump(current, fh)
    return {"rounds": len(summary), "incorporated": len(incorporated),
            "gff3_path": os.path.join(out_dir, "annotation.updated.gff3")}


def run(run: Run, s: Settings, prep: Dict[str, Any], coding: Dict[str, Any]) -> Dict[str, Any]:
    inputs = {"annotation": prep["annotation_models_path"], "isoforms": coding["isoforms_path"],
              "rounds": s.max_rounds, "run_label": s.run_label, "protect": s.protect_supported_models, "additive": s.additive,
              "trust_fl": s.trust_fl_status, "nmd": (s.reject_nmd_targets, s.nmd_stop_to_junction),
              "pct_coding": s.min_percent_coding, "fl_orf": s.min_fl_orf_aa, "min_tpm": s.min_tpm, "min_frac": s.min_isoform_fraction,
              "tss_sites": [run.fdigest(x) for x in s.tss_sites], "polya_sites": [run.fdigest(x) for x in s.polya_sites],
              "site_window": s.site_window, "alt_end_isoforms": s.alt_end_isoforms, "noncoding": s.include_noncoding, "sample_frac": s.min_sample_fraction, "refframe": s.reference_frame_check,
              "code": s.genetic_code, "stage": "annot_update/2", "src": code_digest("pasa_xg.stages.annot_update", "pasa_xg.annot", "pasa_xg.models", "pasa_xg.genome", "pasa_xg.gff")}

    def body(d: str) -> Dict[str, Any]:
        with open(prep["annotation_models_path"], "rb") as fh:
            ann = pickle.load(fh)
        isos = iso_mod.load(coding["isoforms_path"])
        return update(ann, isos, Genome(prep["genome_path"]), s, d,
                      Params(PROTECT_SUPPORTED_MODELS=s.protect_supported_models,
                             ADDITIVE=s.additive,
                             TRUST_FL_STATUS=s.trust_fl_status,
                             REJECT_NMD_TARGETS=s.reject_nmd_targets,
                             NMD_STOP_TO_JUNCTION=s.nmd_stop_to_junction,
                             MIN_PERCENT_PROT_CODING=s.min_percent_coding,
                             MIN_FL_ORF_SIZE=s.min_fl_orf_aa,
                             ALT_END_ISOFORMS=s.alt_end_isoforms,
                             MIN_SAMPLE_FRACTION=s.min_sample_fraction))

    return run.run_stage("annot_update", "05_annot_update", inputs, body)
