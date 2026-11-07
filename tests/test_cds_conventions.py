"""CDS conventions: reference CDS updates (--reference_cds_updates, --additive), NMD candidates keep
their CDS, immunoglobulin / T-cell receptor gene segments are left alone."""

import pandas as pd

from pasa_xg.annot import refframe
from pasa_xg.annot.apply import apply_round
from pasa_xg.annot.checks import Params, cds_extends, nmd_distance
from pasa_xg.annot.noncoding import run_noncoding
from pasa_xg.models import AnnotationSet, Transcript

from test_annot import TABLE, Locus, decide, iso, model


def _upstream_atg_locus(tmp_path):
    """A gene whose reference CDS starts at an internal ATG, 10 codons into the longest ORF: the
    transcript holds an upstream in-frame start codon."""
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 300, 300])
    seq = "".join(L.seq)
    at = ex[0][0] - 1 + 60 + 30          # 0-based: utr5 60, then codon 11 of the ORF
    L.seq = [seq[:at] + "ATG" + seq[at + 3:]]
    L.pos = len(L.seq[0])
    g = L.genome()
    gene = model("G1", ex, g)
    t = gene.transcripts["G1.t1"]
    long_cds = list(t.cds)
    t.cds = [(t.cds[0][0] + 30, t.cds[0][1])] + t.cds[1:]
    assert g.protein(t, TABLE).startswith("M")
    return ex, g, AnnotationSet([gene]), long_cds


def test_cds_extension_test():
    old = Transcript("o", "g", "c", "+", [(1, 300)], cds=[(31, 120)])
    assert not cds_extends(old, old.copy(cds=[(61, 120)]), _G(), TABLE)    # shorter
    assert not cds_extends(old, old.copy(cds=[(32, 121)]), _G(), TABLE)    # reframed
    assert cds_extends(old, old.copy(cds=[(1, 120)]), _G(), TABLE)          # N-terminal extension


class _G:
    """Sequence stub: one residue per reading frame, so proteins compare by length and frame only."""
    def protein(self, t, table):
        n = sum(b - a + 1 for a, b in t.cds)
        return "ACD"[t.cds[0][0] % 3] * (n // 3)


def test_upstream_in_frame_start_extends_the_reference_cds_by_default(tmp_path):
    ex, g, ann, long_cds = _upstream_atg_locus(tmp_path)
    full = iso("full", "L1", ex, g)
    d, decisions, splits = decide(ann, [full], g)
    assert d["full"].action == "protein_update" and d["full"].status == 6
    new, _, _ = apply_round(ann, decisions, splits, "r1", {})
    assert new.genes["G1"].transcripts["G1.t1"].cds == long_cds and len(new.genes["G1"].transcripts) == 1


def test_additive_or_never_adds_the_longer_cds_as_a_new_isoform(tmp_path):
    ex, g, ann, long_cds = _upstream_atg_locus(tmp_path)
    ref_cds = list(ann.genes["G1"].transcripts["G1.t1"].cds)
    for p in (Params(ADDITIVE=True), Params(REFERENCE_CDS_UPDATES="never")):
        d, decisions, splits = decide(ann, [iso("full", "L1", ex, g)], g, p)
        assert d["full"].status == 49 and d["full"].action == "new_isoform"
        new, applied, _ = apply_round(ann, decisions, splits, "r1", {})
        tx = new.genes["G1"].transcripts
        assert tx["G1.t1"].cds == ref_cds                      # the reference CDS is untouched
        assert tx[applied["full"]].exons == ex and tx[applied["full"]].cds == long_cds


def test_reference_start_opening_an_nmd_orf_keeps_it_as_an_nmd_candidate(tmp_path):
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 300, 300, 300])
    g = L.genome()
    ann = AnnotationSet([model("G1", ex, g)])
    shifted = [ex[0], (ex[1][0] + 1, ex[1][1]), ex[2], ex[3]]   # frame shift after exon 1
    i = iso("fs", "L", shifted, g)
    own = i.tx.segments_for(700, 790)
    i.orf, i.tx.cds = None, own
    from pasa_xg.isoforms import Orf
    i.orf = Orf("fs.p1", "complete", 700, 790, 30)
    refframe.check({"fs": i}, ann, g, TABLE, None, Params(), min_aa=10)
    assert i.tx.attrs["pasa_xg_orf"] == "reference_start_nmd" and i.tx.cds
    assert nmd_distance(i.tx) > 50 and i.tx.attrs["transcript_type"] == "nonsense_mediated_decay"
    # the coding rounds reject an NMD target; the noncoding pass adds it with its CDS
    new, d, _ = run_noncoding(ann, {"fs": i}, ["fs"], g, TABLE, Params(), None, require_ends=False)
    t = new.genes["G1"].transcripts[d[0].annotation_id]
    assert t.cds and t.attrs["pasa_xg_coding"] == "nmd" and t.attrs["transcript_type"] == "nonsense_mediated_decay"


def test_ig_tr_gene_segments_are_left_as_annotated(tmp_path):
    from pasa_xg.config import Settings
    from pasa_xg.stages.annot_update import update
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 150, 150, 300])
    g = L.genome()
    ann = AnnotationSet([model("G1", ex, g)])
    ann.genes["G1"].attrs["gene_type"] = "IG_C_gene"
    before = ann.genes["G1"].transcripts["G1.t1"].copy()
    isos = {"skip": iso("skip", "L1", [ex[0], ex[2], ex[3]], g)}
    s = Settings(genome=str(tmp_path / "g.fa"), input_type="transcripts", out_dir=str(tmp_path), reads=["r.fa"])
    update(ann, isos, g, s, str(tmp_path / "upd"), Params())
    d = pd.read_csv(tmp_path / "upd" / "annot_update.decisions.tsv", sep="\t").set_index("isoform_id")
    assert d.loc["skip", "pasa_status"] == 50 and not d.loc["skip", "incorporated"]
    import pickle
    with open(tmp_path / "upd" / "annotation.updated.pkl", "rb") as fh:
        out = pickle.load(fh)
    tx = out.genes["G1"].transcripts
    assert list(tx) == ["G1.t1"] and tx["G1.t1"].exons == before.exons and tx["G1.t1"].cds == before.cds
