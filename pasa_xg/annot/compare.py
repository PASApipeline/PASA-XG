"""One round of PASA's annotation comparison, over LRAA isoforms.

A port of cDNA_annotation_comparer.dbi. PASA's units map as follows:
  - a PASA subcluster (overlapping assemblies of one gene) is an LRAA gene;
  - a full-length (FL) assembly is an isoform with a complete TransDecoder ORF;
    every other isoform takes PASA's EST (non-FL) path;
  - the annotation is fixed during a round; accepted changes are applied afterwards
    (apply.py), and the next round compares against the result.

Each isoform gets a PASA status id (kept for traceability; DESIGN.md §5.2 maps them
to actions) and, when it changes the annotation, the proposed structure.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Set, Tuple

from ..genome import Genome
from ..isoforms import Isoform
from ..models import AnnotationSet, Gene, Span, Transcript
from . import checks as C
from .structure import gene_from_exons, merge_compatible, stitch

log = logging.getLogger("pasa_xg")

STATUS = {
    1: "fl-cdna assembly fails validation tests",
    3: "gene-compatible fl-cdna assembly incorporated",
    4: "gene-compatible fl-cdna assembly alters UTRs",
    5: "fl-cdna assembly matches multiple genes",
    6: "gene-compatible fl-cdna assembly alters protein, passed validation",
    7: "gene-compatible fl-cdna assembly alters protein, fails validation",
    8: "incompatible fl-cdna assembly alignment updates gene structure",
    9: "incompatible fl-cdna assembly provides alternative splicing isoform, passes validation",
    10: "fl-cdna assembly provides a novel gene",
    11: "EST assembly matches multiple genes",
    12: "EST assembly incorporated",
    13: "EST assembly extends UTRs",
    14: "EST assembly alters protein sequence, passes validation",
    15: "EST assembly alters protein sequence, fails validation",
    16: "EST assembly properly stitched into gene structure",
    18: "EST assembly stitched into gene model fails validation",
    19: "EST assembly aligns to intergenic region",
    21: "fl-cdna assembly provides alternate splicing isoform, fails validation",
    22: "fl-cdna assembly provides novel gene which overlaps existing annotations",
    23: "incompatible fl-cdna assembly fails validation",
    24: "FL-cDNA assembly stitched by EST assembly to provide alt splicing isoform",
    25: "EST assembly stitched into gene model requires alternative splicing isoform",
    26: "FL-cDNA spans single gene and allowed to STOMP it",
    27: "EST-assembly stitched into a FL-alignment providing new alt splice isoform",
    28: "EST-assembly stitched into a FL-alignment but fails validation tests",
    29: "FL-cDNA found capable of merging multiple genes",
    30: "FL-cDNA overlaps multiple genes, but found incapable of merging them",
    31: "EST-stitched assembly STOMPS model lacking transcript support",
    32: "EST-stitched gene w/preexisting transcript support STOMPS a new alt splicing variation",
    33: "FL-assembly STOMPS new splice isoform",
    34: "single-exon EST-assembly fails gene compatibility test",
    35: "EST-assembly overlaps multiple genes, but found incapable of merging them",
    36: "EST-assembly found capable of merging multiple genes",
    37: "delayed incorporation of FL-assembly due to successful merging operation",
    38: "delayed incorporation of EST-assembly due to successful merging operation",
    39: "delayed incorporation of FL-assembly due to successful splitting operation",
    40: "FL-cDNAs split single gene into multiple genes",
    41: "FL-cDNAs suggest split single gene into multiple genes, failed automation",
    42: "FL-cDNA aligns as antisense to existing annotated gene",
    43: "EST assembly aligns as antisense to existing annotated gene",
    44: "FL-cDNA provides alt splicing isoform of a novel gene",
    45: "fl-cdna assembly provides isoform for novel gene which overlaps existing annotations",
    100: "below the expression threshold",
}

ACTION = {
    3: "matches_existing", 12: "matches_existing",
    4: "utr_update", 13: "utr_update",
    6: "protein_update", 8: "protein_update", 14: "protein_update", 16: "protein_update",
    9: "new_isoform", 24: "new_isoform", 25: "new_isoform", 27: "new_isoform",
    33: "new_isoform", 44: "new_isoform", 32: "new_isoform",
    29: "merge_genes", 36: "merge_genes",
    40: "split_gene",
    26: "replace_unsupported_model", 31: "replace_unsupported_model",
    10: "novel_gene",
    42: "antisense", 43: "antisense",
    19: "intergenic_noncoding",
    37: "deferred", 38: "deferred", 39: "deferred",
}


def action_of(status: int, valid: bool) -> str:
    a = ACTION.get(status)
    if a in ("antisense", "intergenic_noncoding", "deferred"):
        return a
    if a and valid:
        return a
    return "rejected"


@dataclass
class Cand:
    """One isoform under evaluation (PASA's assembly wrapper)."""
    iso: Isoform
    group: str
    fl: bool
    gene_obj: Optional[Transcript] = None      # the isoform as a coding model
    updated: Optional[Transcript] = None       # proposed structure
    status: Optional[int] = None
    comment: str = ""
    genes: List[str] = field(default_factory=list)
    model: Optional[str] = None                # reference transcript id (or ids for merges)
    valid: bool = False
    alt_splice: bool = False
    novel: bool = False
    novel_gene: Optional[str] = None           # temporary novel gene key

    @property
    def id(self) -> str:
        return self.iso.id

    @property
    def exons(self) -> List[Span]:
        return self.iso.tx.exons

    def set(self, status: int, valid: bool, comment: str = "", **kw) -> None:
        self.status, self.valid = status, valid
        if comment:
            self.comment = (self.comment + "; " + comment).strip("; ") if self.comment else comment
        for k, v in kw.items():
            setattr(self, k, v)


@dataclass
class Decision:
    isoform_id: str
    group: str
    status: int
    action: str
    valid: bool
    genes: List[str]
    model: Optional[str]
    alt_splice: bool
    novel: bool
    novel_gene: Optional[str]
    structure: Optional[Transcript]
    comment: str


class RoundComparer:
    def __init__(self, ann: AnnotationSet, genome: Genome, table: Dict[str, str], p: C.Params,
                 prior_support: Optional[Dict[str, Set[str]]] = None,
                 prior_fl_models: Optional[Set[str]] = None,
                 supported_models: Optional[Set[str]] = None) -> None:
        """`prior_support` / `prior_fl_models`: models built or confirmed by isoforms in
        earlier rounds. PASA forgets these between runs, so a later run may overwrite a
        model an earlier one built from a full-length transcript; rounds remember."""
        self.ann, self.genome, self.table, self.p = ann, genome, table, p
        self.prior_fl_models: Set[str] = set(prior_fl_models or ())
        # models whose structure some isoform supports (see Params.PROTECT_SUPPORTED_MODELS)
        self.protected: Set[str] = set(supported_models or ()) if p.PROTECT_SUPPORTED_MODELS else set()
        if p.ADDITIVE:  # every existing model is kept as is
            self.protected |= {t.id for t in ann.transcripts()}
        self.models_unavail: Set[str] = set(self.prior_fl_models)    # FL-updated models
        self.model_support: Dict[str, Set[str]] = defaultdict(set)   # model -> incorporated isoforms
        for mid, isos in (prior_support or {}).items():
            self.model_support[mid] |= set(isos)
        self.merged_models: Set[str] = set()
        self.novel_by_iso: Dict[str, str] = {}
        self.novel_fl: Dict[str, List[Cand]] = defaultdict(list)     # group -> novel FL cands
        self.validating_fl: Set[str] = set()
        self.demoted: Set[str] = set()                               # FL treated as non-FL
        self.antisense: Set[str] = set()
        self.novel_turned_antisense: Set[str] = set()
        self._novel_counter = 0
        self._partials: Dict[str, Tuple[bool, bool]] = {}

    # ------------------------------------------------------------ helpers
    def partials(self, t: Transcript) -> Tuple[bool, bool]:
        if t.id not in self._partials:
            self._partials[t.id] = C.partial_flags(t, self.genome, self.table)
        return self._partials[t.id]

    def gene(self, exons, strand, contig, partials=(False, False), tid="tmp") -> Transcript:
        return gene_from_exons(exons, strand, contig, self.genome, self.table, tid=tid,
                               allow_5partial=partials[0], allow_3partial=partials[1])

    def compare(self, old: Transcript, new: Transcript, fl_alt: bool = False,
                new_partials: Tuple[bool, bool] = (False, False)) -> Tuple[bool, str]:
        return C.compare_proteins(old, new, self.genome, self.table, self.p, fl_alt, *new_partials)

    def mapped_genes(self, cands: List[Cand], fli: bool, strand: Optional[str]) -> List[str]:
        """PASA's get_potential_overlapping_annotations for a set of isoforms."""
        if not cands:
            return []
        contig = cands[0].iso.tx.contig
        lend = min(c.iso.tx.start for c in cands)
        rend = max(c.iso.tx.end for c in cands)
        segs = [e for c in cands for e in c.exons]
        out = []
        for g in self.ann.overlapping(contig, lend, rend):
            if strand is not None and g.strand != strand:
                continue
            if not any(C.span_overlap(e, s) for t in g.transcripts.values() for e in t.exons for s in segs):
                continue
            common = C.span_overlap((lend, rend), (g.start, g.end))
            if (not fli or 100.0 * common / (g.end - g.start + 1) >= self.p.MIN_PERCENT_OVERLAP
                    or 100.0 * common / (rend - lend + 1) >= self.p.MIN_PERCENT_OVERLAP):
                out.append(g.id)
        return out

    def models(self, gene_id: str) -> List[Transcript]:
        return [t for t in self.ann.genes[gene_id].transcripts.values() if t.cds]

    def new_novel_gene(self) -> str:
        self._novel_counter += 1
        return f"novel_{self._novel_counter}"

    # ------------------------------------------------------------ FL path
    def validate_fl(self, cands: List[Cand]) -> Tuple[List[Cand], List[Cand], List[Cand]]:
        """Returns (validating, still FL, demoted to non-FL)."""
        validating, fl, demoted = [], [], []
        for c in cands:
            t = c.iso.tx
            if not t.cds:
                demoted.append(c)
                continue
            fl.append(c)
            g = t.copy(id=c.id)
            c.gene_obj = c.updated = g
            errs = C.validate_gene(g, self.genome, self.table, self.p)
            pep_len = g.cds_length() // 3
            pct_coding = 100.0 * pep_len * 3 / g.length()
            if pct_coding < self.p.MIN_PERCENT_PROT_CODING:
                if not (self.p.MIN_FL_ORF_SIZE and pep_len > self.p.MIN_FL_ORF_SIZE):
                    errs.append(f"protein coding only {pct_coding:.0f}% of the transcript")
            if not errs:
                ok, msg = C.check_utr_exons(g, self.p)
                if not ok:
                    errs.append(msg)
            if errs:
                c.set(1, False, "; ".join(errs))
            else:
                self.validating_fl.add(c.id)
                validating.append(c)
        return validating, fl, demoted

    def analyze_fl(self, cands: List[Cand]) -> List[Cand]:
        """FL isoforms of one group. Returns those to treat as non-FL."""
        validating, fl, demoted = self.validate_fl(cands)
        if validating:
            validating.sort(key=lambda c: (c.gene_obj.num_exons, c.gene_obj.cds_length()), reverse=True)
            demoted += self.improve_with_fl(validating)
        novel = any(c.novel for c in validating)
        genes_text = validating[0].genes if validating else []
        if genes_text and not novel:
            for c in cands:
                if not c.genes and c not in demoted:
                    c.genes = list(genes_text)
        else:
            usable = [c for c in fl if c not in demoted]
            strand = cands[0].iso.tx.strand
            hits = self.mapped_genes(usable, fli=False, strand=strand)
            if hits:
                for c in cands:
                    if not c.genes and c not in demoted:
                        c.genes = list(hits)
                    if c.status in (10, 44):
                        c.set(22 if c.status == 10 else 45, False,
                              "novel gene overlaps an existing annotation")
                        self.novel_fl[c.group] = [x for x in self.novel_fl[c.group] if x is not c]
        for c in demoted:
            self.demoted.add(c.id)
            c.fl = False
            c.status, c.valid, c.updated = None, False, None
        return demoted

    def improve_with_fl(self, valid: List[Cand]) -> List[Cand]:
        demote: List[Cand] = []
        strand = valid[0].iso.tx.strand
        genes = self.mapped_genes(valid, fli=True, strand=strand)
        if len(genes) > 1:
            for c in valid:
                c.genes = list(genes)
                c.set(5, False, "matches multiple genes")
            return demote
        if not genes:
            key = self.new_novel_gene()
            for i, c in enumerate(valid):
                c.genes = [key]
                c.novel, c.novel_gene = True, key
                self.novel_by_iso[c.id] = key
                self.novel_fl[c.group].append(c)
                if i == 0:
                    c.set(10, True, "no overlapping gene: novel")
                else:
                    c.set(44, True, "alternative isoform of a novel gene", alt_splice=True)
            return demote
        gid = genes[0]
        for c in valid:
            c.genes = [gid]
        models = self.models(gid)
        by_id = {m.id: m for m in models}
        model_to_cand: Dict[str, str] = {}
        cand_to_model: Dict[str, str] = {}
        demoted_ids: Set[str] = set()
        # associate compatible models and FL isoforms
        for m in models:
            p5, p3 = self.partials(m)
            for c in valid:
                if c.id in demoted_ids:
                    continue
                if m.id in model_to_cand:
                    break
                if c.id in cand_to_model:
                    continue
                merged = merge_compatible(m.exons, c.exons)
                if merged is None:
                    continue
                if (c.iso.tx.start > m.start or c.iso.tx.end < m.end) and not self.p.TRUST_FL_STATUS:
                    demote.append(c)
                    demoted_ids.add(c.id)
                    c.comment = "compatible FL isoform is extended by the gene model; treated as non-FL"
                    continue
                model_to_cand[m.id], cand_to_model[c.id] = c.id, m.id
                mg = self.gene(merged, m.strand, m.contig, (p5, p3), tid=c.id)
                ok, _ = self.compare(m, mg, new_partials=(p5, p3))
                if ok:
                    c.updated = mg
        trying = [c for c in valid if c.id not in demoted_ids]
        # models built from full-length isoforms in earlier rounds are not replaced
        updated: Set[str] = {m.id for m in models if m.id in self.prior_fl_models}
        protected = {m.id for m in models if m.id in self.protected}
        # compatible updates
        for c in trying:
            mid = cand_to_model.get(c.id)
            if not mid:
                continue
            m, new = by_id[mid], c.updated
            if C.same_exons(new.exons, m.exons) and new.cds == m.cds:
                c.set(3, True, "identical to the gene model", model=mid)
            elif new.cds == m.cds:
                c.set(4, True, "CDS same, UTRs differ", model=mid)
            else:
                ok, msg = self.compare(m, new, new_partials=self.partials(m))
                c.set(6 if ok else 7, ok, msg, model=mid)
            if c.valid:
                updated.add(mid)
                self.model_support[mid].add(c.id)
        # incompatible updates of models not yet updated (nor protected)
        available = [m for m in models if m.id not in updated and m.id not in protected]
        if available:
            for c in trying:
                if c.status:
                    continue
                for m in available:
                    if m.id in updated:
                        continue
                    ok, msg = self.compare(m, c.gene_obj)
                    c.comment, c.model = msg, m.id
                    if ok:
                        c.set(8, True, "replaces the incompatible gene model", updated=c.gene_obj)
                        updated.add(m.id)
                        self.model_support[m.id].add(c.id)
                        break
            if self.p.STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE:
                for c in trying:
                    if c.status or c.gene_obj.num_exons == 1:
                        continue
                    for m in available:
                        if m.id in updated:
                            continue
                        ok, msg = self.stomp_ok(m, c.gene_obj)
                        if ok:
                            c.set(26, True, msg, model=m.id, updated=c.gene_obj)
                            updated.add(m.id)
                            self.model_support[m.id].add(c.id)
                            break
        # alternative isoforms, once every model is accounted for
        if not [m for m in models if m.id not in updated and m.id not in protected]:
            for c in trying:
                if c.status:
                    continue
                ok, msg = False, ""
                for m in models:
                    ok, msg = self.compare(m, c.gene_obj, fl_alt=True)
                    if ok:
                        c.set(9, True, "new alternative isoform; " + msg, model=m.id,
                              alt_splice=True, updated=c.gene_obj)
                        break
                if not ok and self.p.STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE:
                    for m in models:
                        ok, smsg = self.stomp_ok(m, c.gene_obj)
                        if ok:
                            c.set(33, True, smsg, model=m.id, alt_splice=True, updated=c.gene_obj)
                            break
                if not ok:
                    c.set(21, False, msg, alt_splice=True)
        for c in trying:
            if not c.status:
                c.set(23, False, "could not be incorporated")
        self.models_unavail |= updated
        return demote

    def stomp_ok(self, model: Transcript, new: Transcript) -> Tuple[bool, str]:
        a, b = C.cds_span(model), C.cds_span(new)
        pct = 100.0 * C.span_overlap(a, b) / (a[1] - a[0] + 1)
        ok = pct >= self.p.MIN_PERCENT_OVERLAP_GENE_REPLACE
        return ok, f"stomp: {pct:.0f}% overlap of {model.id}"

    # ------------------------------------------------------------ non-FL path
    def analyze_nonfl(self, cands: List[Cand]) -> None:
        group = cands[0].group
        if self.novel_fl.get(group):
            self.stitch_into_novel(cands, self.novel_fl[group])
            return
        strand = cands[0].iso.tx.strand
        genes = self.mapped_genes(cands, fli=False, strand=strand)
        for c in cands:
            c.genes = list(genes)
        if len(genes) > 1:
            for c in cands:
                c.set(11, False, "matches multiple genes")
            return
        if not genes:
            for c in cands:
                c.set(19, False, "no overlapping gene")
            return
        self.incorporate_nonfl(cands, genes[0])

    def incorporate_nonfl(self, cands: List[Cand], gid: str) -> None:
        models = self.models(gid)
        merged_into: Dict[str, str] = {}
        original = {m.id: m for m in models}
        current_exons: Dict[str, List[Span]] = {m.id: list(m.exons) for m in models}
        current_gene: Dict[str, Transcript] = {m.id: m for m in models}
        available = [m.id for m in models if m.id not in self.models_unavail]
        fl_updated = [m.id for m in models if m.id in self.models_unavail]
        for mid in available:
            m = original[mid]
            parts = self.partials(m)
            for c in cands:
                if c.id in merged_into:
                    continue
                merged = merge_compatible(current_exons[mid], c.exons)
                if merged is None:
                    continue
                new = self.gene(merged, m.strand, m.contig, parts, tid=mid)
                incorporated = False
                if C.same_exons(new.exons, m.exons) and new.cds == m.cds:
                    status, msg, incorporated = 12, "identical to the gene model", True
                elif new.cds == m.cds:
                    status, msg, incorporated = 13, "CDS same, UTRs differ", True
                else:
                    ok, msg = self.compare(m, new, new_partials=parts)
                    if ok:
                        ok, umsg = C.check_utr_exons(new, self.p)
                        msg = msg if ok else umsg
                    status, incorporated = (14, True) if ok else (15, False)
                if incorporated:
                    current_gene[mid], current_exons[mid] = new, merged
                    merged_into[c.id] = mid
                    self.model_support[mid].add(c.id)
                    c.set(status, True, msg, model=mid)
                else:
                    c.set(status, False, msg, model=mid)
        for c in cands:
            if c.id not in merged_into and c.iso.tx.num_exons == 1:
                c.set(34, False, "single-exon isoform incompatible with the gene models")
        # stitching
        for c in cands:
            if c.id in merged_into or c.status:
                continue
            valid, msg, templates = False, "", []
            have_fl_stitched, needs_isoform = False, False
            last_example = None
            for mid in available + fl_updated:
                m = original[mid]
                parts = self.partials(m)
                st_exons = stitch(current_exons[mid], c.exons)
                st = self.gene(st_exons, m.strand, m.contig, parts, tid=mid)
                last_example = (mid, st, st_exons)
                valid, msg = C.check_false_utrs((c.iso.tx.start, c.iso.tx.end), m, st)
                if valid:
                    valid, msg = self.compare(m, st, new_partials=parts)
                if valid:
                    valid, umsg = C.check_utr_exons(st, self.p)
                    msg = msg if valid else umsg
                if valid:
                    templates.append((mid, st))
                    if mid in self.models_unavail:
                        have_fl_stitched = True
                if valid and not have_fl_stitched and mid in self.protected:
                    needs_isoform = True
                    msg += "; the model's own structure is supported, so this is a new isoform"
                elif valid and not have_fl_stitched:
                    supported = [x for x in cands if merged_into.get(x.id) == mid]
                    final = st_exons
                    for x in supported:
                        final = merge_compatible(final, x.exons) if final is not None else None
                    if supported and final is None:
                        needs_isoform = True
                        msg += "; stitched model incompatible with isoforms already incorporated"
                    else:
                        if final is not None and final != st_exons:
                            st = self.gene(final, m.strand, m.contig, parts, tid=mid)
                        merged_into[c.id] = mid
                        current_gene[mid], current_exons[mid] = st, final or st_exons
                        self.model_support[mid].add(c.id)
                if valid and not have_fl_stitched:
                    break
            if valid and not have_fl_stitched and not needs_isoform:
                c.set(16, True, msg, model=merged_into[c.id])
            elif valid and (needs_isoform or (have_fl_stitched and templates)):
                if have_fl_stitched:
                    mid, st = max(templates, key=lambda x: x[1].cds_length())
                    c.set(24, True, msg, model=mid, updated=st, alt_splice=True)
                else:
                    mid, st, _ = last_example
                    c.set(25, True, msg, model=mid, updated=st, alt_splice=True)
            else:
                c.set(18, False, (msg + "; " if msg else "") + "stitched model fails validation",
                      model=last_example[0] if last_example else None,
                      updated=last_example[1] if last_example else None)
        if self.p.STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE:
            self.stomp_nonfl(cands, original, current_exons, current_gene, merged_into,
                             available, fl_updated)
        # valid single-model updates carry the model's final structure
        for c in cands:
            if c.valid and not c.alt_splice and c.model in current_gene:
                c.updated = current_gene[c.model]

    def stomp_nonfl(self, cands, original, current_exons, current_gene, merged_into,
                    available, fl_updated) -> None:
        for c in cands:
            if c.status == 34 or c.id in merged_into:
                continue
            unsupported = [m for m in available if not self.model_support.get(m)]
            supported = [m for m in available if self.model_support.get(m)]
            for mid in unsupported + fl_updated + supported:
                m = current_gene[mid]
                parts = self.partials(original[mid])
                st_exons = stitch(current_exons[mid], c.exons)
                st = self.gene(st_exons, m.strand, m.contig, parts, tid=mid)
                ok, _ = C.check_false_utrs((c.iso.tx.start, c.iso.tx.end), m, st)
                ok = ok and C.check_utr_exons(st, self.p)[0]
                if not ok:
                    continue
                ok, msg = self.stomp_ok(m, st)
                if ok:
                    merged_into[c.id] = mid
                    has_support = bool(self.model_support.get(mid)) or mid in fl_updated or self.p.ADDITIVE
                    if not has_support:
                        c.set(31, True, msg, model=mid, updated=st)
                        current_gene[mid], current_exons[mid] = st, st_exons
                        self.model_support[mid].add(c.id)
                    else:
                        c.set(32, True, msg, model=mid, updated=st, alt_splice=True)
                    break

    def stitch_into_novel(self, cands: List[Cand], novel: List[Cand]) -> None:
        for c in cands:
            for fl in novel:
                st_exons = stitch(fl.gene_obj.exons, c.exons)
                st = self.gene(st_exons, fl.gene_obj.strand, fl.gene_obj.contig)
                ok, msg = C.check_utr_exons(st, self.p)
                if ok:
                    ok, msg = self.compare(fl.gene_obj, st)
                c.genes, c.novel, c.novel_gene = [fl.novel_gene], True, fl.novel_gene
                c.set(27 if ok else 28, ok, msg, alt_splice=True, updated=st)
                if ok:
                    break

    # ------------------------------------------------------------ multi-gene, antisense, split
    def analyze_multi_gene(self, cands: List[Cand], is_fl: bool) -> List[Cand]:
        """Returns isoforms whose incorporation must wait (merged genes' support)."""
        deferred: List[Cand] = []
        for c in cands:
            if self.p.ADDITIVE:
                c.set(30 if is_fl else 35, False, "additive update: existing genes are not merged")
                continue
            genes = [self.ann.genes[g] for g in c.genes]
            if any(len([t for t in g.transcripts.values()]) > 1 for g in genes):
                c.set(30 if is_fl else 35, False, "a matched gene already has alternative isoforms")
                continue
            models = [next(iter(g.transcripts.values())) for g in genes]
            if any(not m.cds for m in models):
                c.set(30 if is_fl else 35, False, "a matched gene is noncoding")
                continue
            if any(m.id in self.models_unavail or (self.model_support.get(m.id) and any(
                    x in self.validating_fl for x in self.model_support[m.id])) for m in models):
                c.set(30 if is_fl else 35, False, "a matched gene already has full-length support")
                continue
            if len({g.strand for g in genes}) != 1:
                c.set(30 if is_fl else 35, False, "matched genes are on opposite strands")
                continue
            if is_fl:
                update = c.gene_obj
                ok, msg = True, ""
            else:
                st_exons = c.exons
                for m in models:
                    st_exons = stitch(m.exons, st_exons)
                update = self.gene(st_exons, models[0].strand, models[0].contig)
                ok, msg = True, ""
                for m in models:
                    ok, msg = self.compare(m, update)
                    if not ok:
                        break
                if ok:
                    combined = Transcript("combined", "x", models[0].contig, models[0].strand,
                                          _union_exons([e for m in models for e in m.exons]),
                                          cds=[s for m in models for s in m.cds])
                    ok, msg = C.check_false_utrs((c.iso.tx.start, c.iso.tx.end), combined, update)
                if ok:
                    ok, msg = C.check_utr_exons(update, self.p)
            c.updated = update
            if not ok:
                c.set(30 if is_fl else 35, False, msg)
                continue
            ul, ur = C.cds_span(update)
            fail = ""
            for m in models:
                if m.id in self.merged_models:
                    fail = f"{m.id} already merged"
                    break
                ml, mr = C.cds_span(m)
                cov = 100.0 * C.span_overlap((ul, ur), (ml, mr)) / (mr - ml + 1)
                if cov < self.p.MIN_PERCENT_OVERLAP_GENE_REPLACE:
                    fail = f"{m.id} covered only {cov:.0f}%"
                    break
            if fail:
                c.set(30 if is_fl else 35, False, fail)
                continue
            for m in models:
                self.merged_models.add(m.id)
                for iso_id in self.model_support.pop(m.id, set()):
                    deferred.append(iso_id)
            c.set(29 if is_fl else 36, True, "merges " + " ".join(g.id for g in genes),
                  model=" ".join(m.id for m in models))
        return deferred

    def classify_antisense(self, cands: Iterable[Cand]) -> None:
        for c in cands:
            if c.status not in (1, 10, 44, 19, 22, 27, 34):
                continue
            flip = "-" if c.iso.tx.strand == "+" else "+"
            genes = self.mapped_genes([c], fli=False, strand=flip)
            if genes:
                if c.novel_gene and c.status in (10, 22):
                    self.novel_turned_antisense.add(c.novel_gene)
                c.genes = genes
                c.set(42 if c.fl else 43, False, "antisense to " + " ".join(genes))
                self.antisense.add(c.id)

    def analyze_splitting(self, cands: List[Cand]) -> Tuple[List[Tuple[str, List[Cand]]], List[str]]:
        """Genes that FL isoforms from two or more LRAA genes split. Returns the
        approved splits (gene id, isoforms) and isoform ids to defer."""
        by_gene: Dict[str, Dict[str, List[Cand]]] = defaultdict(lambda: defaultdict(list))
        for c in cands:
            if (c.id in self.validating_fl and c.fl and c.id not in self.antisense
                    and len(c.genes) == 1 and c.genes[0] in self.ann.genes and not c.novel):
                by_gene[c.genes[0]][c.group].append(c)
        splits, deferred = [], []
        if self.p.ADDITIVE:  # existing genes are not split
            return splits, deferred
        for gid, groups in by_gene.items():
            if len(groups) < 2:
                continue
            g = self.ann.genes[gid]
            if len(g.transcripts) != 1:
                continue
            m = next(iter(g.transcripts.values()))
            if not m.cds or m.id in self.models_unavail or m.id in self.merged_models:
                continue
            ml, mr = C.cds_span(m)
            chosen = []
            for grp, cs in groups.items():
                cs = [c for c in cs if C.span_overlap(C.cds_span(c.gene_obj), (ml, mr))]
                if cs:
                    chosen.append(max(cs, key=lambda c: c.iso.tx.length()))
            if len(chosen) < 2 or any(c.iso.tx.strand != m.strand for c in chosen):
                continue
            lo = min(c.iso.tx.start for c in chosen)
            hi = max(c.iso.tx.end for c in chosen)
            cov = 100.0 * C.span_overlap((lo, hi), (ml, mr)) / (mr - ml + 1)
            if cov >= self.p.MIN_PERCENT_OVERLAP_GENE_REPLACE:
                deferred += list(self.model_support.pop(m.id, set()))
                for c in chosen:
                    c.set(40, True, f"splits {gid}", model=m.id, updated=c.gene_obj, alt_splice=False)
                splits.append((gid, chosen))
            else:
                for c in chosen:
                    c.set(41, False, f"suggests splitting {gid}, but covers only {cov:.0f}% of its CDS")
        return splits, deferred


def _union_exons(spans: List[Span]) -> List[Span]:
    out: List[List[int]] = []
    for a, b in sorted(spans):
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def supported_models(ann: AnnotationSet, isoforms: Iterable[Isoform]) -> Set[str]:
    """Models whose structure an isoform supports: a multi-exon model with a compatible
    isoform sharing at least one of its introns; a single-exon model with a compatible
    single-exon isoform."""
    by_contig: Dict[Tuple[str, str], List[Transcript]] = defaultdict(list)
    for i in isoforms:
        by_contig[(i.tx.contig, i.tx.strand)].append(i.tx)
    for v in by_contig.values():
        v.sort(key=lambda t: t.start)
    out: Set[str] = set()
    for g in ann:
        txs = by_contig.get((g.contig, g.strand), [])
        for m in g.transcripts.values():
            m_introns = set(m.introns())
            for t in txs:
                if t.start > m.end:
                    break
                if t.end < m.start:
                    continue
                if m_introns and not (m_introns & set(t.introns())):
                    continue
                if not m_introns and t.num_exons > 1:
                    continue
                if merge_compatible(m.exons, t.exons) is not None:
                    out.add(m.id)
                    break
    return out


def run_round(ann: AnnotationSet, isoforms: List[Isoform], genome: Genome, table: Dict[str, str],
              p: C.Params, prior_support: Optional[Dict[str, Set[str]]] = None,
              prior_fl_models: Optional[Set[str]] = None,
              all_isoforms: Optional[Iterable[Isoform]] = None,
              ) -> Tuple[List[Decision], List[Tuple[str, List[Cand]]], RoundComparer]:
    """Compare `isoforms` with `ann`; returns decisions, approved splits and the comparer.
    Model support is judged against `all_isoforms` (default: `isoforms`)."""
    supported = supported_models(ann, all_isoforms if all_isoforms is not None else isoforms) \
        if p.PROTECT_SUPPORTED_MODELS else set()
    rc = RoundComparer(ann, genome, table, p, prior_support, prior_fl_models, supported)
    groups: Dict[str, List[Cand]] = defaultdict(list)
    for iso in isoforms:
        groups[iso.tx.gene_id].append(Cand(iso, iso.tx.gene_id, fl=iso.full_length))
    order = sorted(groups, key=lambda g: (groups[g][0].iso.tx.contig, min(c.iso.tx.start for c in groups[g])))
    all_cands = [c for g in order for c in groups[g]]
    # FL mode, then non-FL mode (PASA processes each contig's clusters this way)
    nonfl: Dict[str, List[Cand]] = {g: [c for c in groups[g] if not c.fl] for g in order}
    for g in order:
        fl = [c for c in groups[g] if c.fl]
        if fl:
            nonfl[g] += rc.analyze_fl(fl)
    # the most complete candidate gets the first chance to merge
    most_complete = lambda c: (-c.iso.tx.num_exons, -(c.gene_obj or c.iso.tx).cds_length(), c.id)
    multi_fl = sorted((c for c in all_cands if c.status == 5), key=most_complete)
    for g in order:
        if nonfl[g]:
            rc.analyze_nonfl(nonfl[g])
    multi_nonfl = sorted((c for c in all_cands if c.status == 11), key=most_complete)
    deferred_ids = rc.analyze_multi_gene(multi_fl, True) + rc.analyze_multi_gene(multi_nonfl, False)
    rc.classify_antisense(all_cands)
    for c in all_cands:  # alternative isoforms of novel genes that turned out antisense
        if c.novel_gene in rc.novel_turned_antisense and c.status in (27, 44):
            c.set(42 if c.fl else 43, False, "isoform of a novel gene found antisense")
    splits, split_deferred = rc.analyze_splitting(all_cands)
    by_id = {c.id: c for c in all_cands}
    for iso_id in deferred_ids:
        c = by_id.get(iso_id)
        if c and c.valid:
            c.set(37 if c.fl else 38, False, "deferred: its gene was merged")
    for iso_id in split_deferred:
        c = by_id.get(iso_id)
        if c and c.valid and c.status != 40:
            c.set(39, False, "deferred: its gene was split")
    decisions = []
    for c in all_cands:
        if c.status is None:
            c.set(23 if c.fl else 18, False, "not classified")
        decisions.append(Decision(c.id, c.group, c.status, action_of(c.status, c.valid), c.valid,
                                  c.genes, c.model, c.alt_splice, c.novel, c.novel_gene,
                                  c.updated if c.valid else None, c.comment))
    return decisions, splits, rc
