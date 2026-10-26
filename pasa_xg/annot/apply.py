"""Apply one round's accepted decisions to the annotation (PASA's
cDNA_annotation_updater): model updates, splits, merges, novel genes and new
alternative isoforms, in that order."""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

from ..models import AnnotationSet, Gene, Transcript
from .compare import Decision

SINGLE_MODEL_UPDATES = {"utr_update", "protein_update", "replace_unsupported_model"}


class Namer:
    """New ids: <gene_id>.pxg<N> for isoforms of known genes, PXG_gene_<N> for
    new genes (their transcripts PXG_gene_<N>.t<M>)."""

    def __init__(self, ann: AnnotationSet) -> None:
        self.next_gene = 1 + max([int(m.group(1)) for g in ann.genes
                                  if (m := re.fullmatch(r"PXG_gene_(\d+)", g))] or [0])
        self.next_iso: Dict[str, int] = defaultdict(lambda: 1)
        for g in ann:
            for tid in g.transcripts:
                m = re.fullmatch(re.escape(g.id) + r"\.pxg(\d+)", tid)
                if m:
                    self.next_iso[g.id] = max(self.next_iso[g.id], int(m.group(1)) + 1)

    def gene(self) -> str:
        gid = f"PXG_gene_{self.next_gene}"
        self.next_gene += 1
        return gid

    def isoform(self, gene_id: str) -> str:
        n = self.next_iso[gene_id]
        self.next_iso[gene_id] += 1
        return f"{gene_id}.pxg{n}"


def _gene_of_model(ann: AnnotationSet) -> Dict[str, str]:
    return {t.id: g.id for g in ann for t in g.transcripts.values()}


def _same_structure(a: Transcript, b: Transcript) -> bool:
    return a.exons == b.exons and a.cds == b.cds


def apply_round(ann: AnnotationSet, decisions: List[Decision], splits, label: str,
                expression: Dict[str, Dict[str, float]], run_label: str = ""
                ) -> Tuple[AnnotationSet, Dict[str, str], Set[str]]:
    """Returns (new annotation, isoform id -> annotation transcript id it produced or
    matched, ids of genes changed or created). Changed models record the latest update
    (pasa_xg_update/round/isoform) and append it to pasa_xg_history, which carries over
    when a later run (`run_label`) updates the output again."""
    new = ann.copy(label=label)
    namer = Namer(new)
    gene_of = _gene_of_model(new)
    applied: Dict[str, str] = {}
    changed: Set[str] = set()
    gone: Set[str] = set()

    def tag(t: Transcript, d: Decision, kind: str) -> None:
        e = expression.get(d.isoform_id, {})
        event = f"{run_label + '/' if run_label else ''}{label}:{kind}:{d.isoform_id}"
        prev = t.attrs.get("pasa_xg_history")
        t.attrs = {**t.attrs, "pasa_xg_update": kind, "pasa_xg_round": label,
                   "pasa_xg_isoform": d.isoform_id,
                   **({"pasa_xg_run": run_label} if run_label else {}),
                   "pasa_xg_history": f"{prev},{event}" if prev else event,
                   **({"TPM": f"{e['TPM']:.3f}"} if "TPM" in e else {})}

    # model updates (several isoforms may build one model's final structure)
    final: Dict[str, Decision] = {}
    for d in decisions:
        if d.valid and d.action in SINGLE_MODEL_UPDATES and d.model in gene_of and d.structure:
            final[d.model] = d
        elif d.valid and d.action == "matches_existing" and d.model in gene_of:
            applied[d.isoform_id] = d.model
    for mid, d in final.items():
        g = new.genes[gene_of[mid]]
        old = g.transcripts[mid]
        t = d.structure.copy(id=mid, gene_id=g.id, attrs=dict(old.attrs))
        tag(t, d, d.action)
        g.transcripts[mid] = t
        changed.add(g.id)
    for d in decisions:
        if d.valid and d.action in SINGLE_MODEL_UPDATES and d.model in final:
            applied[d.isoform_id] = d.model

    # splits: the first product keeps the gene and model ids
    for gid, cands in splits:
        if gid not in new.genes:
            continue
        g = new.genes[gid]
        mid = next(iter(g.transcripts))
        for i, c in enumerate(sorted(cands, key=lambda c: c.gene_obj.start)):
            d = next(x for x in decisions if x.isoform_id == c.id)
            if i == 0:
                t = c.gene_obj.copy(id=mid, gene_id=gid)
                tag(t, d, "split_gene")
                g.transcripts = {mid: t}
                changed.add(gid)
            else:
                ng = namer.gene()
                t = c.gene_obj.copy(id=f"{ng}.t1", gene_id=ng)
                tag(t, d, "split_gene")
                new.add_gene(Gene(ng, t.contig, t.strand, {ng + ".t1": t}, {"split_from": gid}))
                changed.add(ng)
            applied[c.id] = t.id

    # merges: one gene replaces several; it keeps the leftmost gene's ids
    for d in decisions:
        if not (d.valid and d.action == "merge_genes" and d.structure):
            continue
        genes = [new.genes[gid] for gid in d.genes if gid in new.genes]
        if len(genes) != len(d.genes) or any(g.id in gone for g in genes):
            continue
        genes.sort(key=lambda g: g.start)
        keep = genes[0]
        mid = next(iter(keep.transcripts))
        t = d.structure.copy(id=mid, gene_id=keep.id)
        tag(t, d, "merge_genes")
        for g in genes[1:]:
            new.remove_gene(g.id)
            gone.add(g.id)
        keep.transcripts = {mid: t}
        keep.attrs = {**keep.attrs, "merged_from": ",".join(g.id for g in genes)}
        changed.add(keep.id)
        applied[d.isoform_id] = mid

    # novel genes and their isoforms
    novel_ids: Dict[str, str] = {}
    for d in decisions:
        if d.valid and d.action == "novel_gene" and d.structure and d.novel_gene:
            ng = namer.gene()
            novel_ids[d.novel_gene] = ng
            t = d.structure.copy(id=f"{ng}.t1", gene_id=ng)
            tag(t, d, "novel_gene")
            new.add_gene(Gene(ng, t.contig, t.strand, {t.id: t}))
            changed.add(ng)
            applied[d.isoform_id] = t.id

    # new alternative isoforms
    for d in decisions:
        if not (d.valid and d.action == "new_isoform" and d.structure):
            continue
        if d.novel:
            gid = novel_ids.get(d.novel_gene)
            if gid is None:
                continue
            tid = f"{gid}.t{len(new.genes[gid].transcripts) + 1}"
        else:
            gid = gene_of.get(d.model)
            if gid is None or gid not in new.genes or gid in gone:
                continue
            tid = None
        g = new.genes[gid]
        cand = d.structure.copy(gene_id=gid)
        dup = next((t for t in g.transcripts.values() if _same_structure(t, cand)), None)
        if dup is not None:
            applied[d.isoform_id] = dup.id
            continue
        if cand.strand != g.strand or cand.contig != g.contig:
            continue
        cand.id = tid or namer.isoform(gid)
        tag(cand, d, "new_isoform")
        g.add(cand)
        changed.add(gid)
        applied[d.isoform_id] = cand.id
    return new, applied, changed
