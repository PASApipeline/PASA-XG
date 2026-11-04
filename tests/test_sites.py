"""TSS / polyA-aware update (docs/TSS_POLYA_DESIGN.md)."""

from pasa_xg.annot.apply import apply_round
from pasa_xg.annot.checks import Params
from pasa_xg.annot.compare import run_round
from pasa_xg.annot.sites import POLYA, TSS, SiteIndex, add_isoform_ends
from pasa_xg.models import AnnotationSet

from test_annot import TABLE, Locus, iso, model


def decide(ann, isos, genome, sites=None, p=None):
    decisions, splits, rc = run_round(ann, isos, genome, TABLE, p or Params(), sites=sites)
    return {d.isoform_id: d for d in decisions}, decisions, splits, rc


def sites_for(isos, known=()):
    idx = SiteIndex(50)
    for contig, strand, kind, pos in known:
        idx.add(contig, strand, kind, pos, pos, "known:test")
    add_isoform_ends(idx, isos)
    return idx


def test_site_index_window_and_contig_names():
    idx = SiteIndex(50)
    idx.add("22", "+", TSS, 1000, 1010, "known:cage")
    assert idx.near("chr22", "+", TSS, 1060) == "known:cage"    # 50 bp past the peak's end
    assert idx.near("chr22", "+", TSS, 1061) is None
    assert idx.near("chr22", "-", TSS, 1005) is None             # other strand
    assert idx.near("chr22", "+", POLYA, 1005) is None           # other kind
    assert idx.near("chr22", "+", TSS, 1005, sources=("lraa",)) is None


def _locus(tmp_path):
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([400, 300, 300], utr5=250)    # the CDS starts 250 bp into exon 1
    g = L.genome()
    short5 = [(ex[0][0] + 150, ex[0][1])] + ex[1:]   # model with a 150 bp shorter 5' UTR
    return ex, short5, g


def test_unsupported_end_does_not_move_an_anchored_model_end(tmp_path):
    ex, short5, g = _locus(tmp_path)
    ann = AnnotationSet([model("G1", short5, g)])
    longer = iso("longer", "L1", ex, g)
    d, *_ = decide(ann, [longer], g)                              # no sites: PASA extends
    assert d["longer"].action == "utr_update"
    anchored = sites_for([longer], known=[("c", "+", TSS, short5[0][0])])
    d, *_ = decide(ann, [longer], g, anchored)
    assert d["longer"].action == "matches_existing"


def test_supported_end_extends_an_unanchored_model_end(tmp_path):
    ex, short5, g = _locus(tmp_path)
    ann = AnnotationSet([model("G1", short5, g)])
    longer = iso("longer", "L1", ex, g)
    longer.tss = True                                              # LRAA called its TSS
    d, decisions, splits, rc = decide(ann, [longer], g, sites_for([longer]))
    assert d["longer"].action == "utr_update"
    new, _, _ = apply_round(ann, decisions, splits, "r1", {}, alt_sites=rc.alt_sites)
    assert new.genes["G1"].transcripts["G1.t1"].exons == ex


def test_supported_alt_end_at_an_anchored_model_end_is_recorded(tmp_path):
    ex, short5, g = _locus(tmp_path)
    ann = AnnotationSet([model("G1", short5, g)])
    longer = iso("longer", "L1", ex, g)
    longer.tss = True
    idx = sites_for([longer], known=[("c", "+", TSS, short5[0][0])])
    d, decisions, splits, rc = decide(ann, [longer], g, idx)
    assert d["longer"].action == "matches_existing"
    new, _, _ = apply_round(ann, decisions, splits, "r1", {}, alt_sites=rc.alt_sites)
    t = new.genes["G1"].transcripts["G1.t1"]
    assert t.exons == short5 and t.attrs["alt_tss_sites"] == f"{ex[0][0]}:lraa"


def test_same_protein_end_variant_is_not_a_new_isoform(tmp_path):
    ex, _, g = _locus(tmp_path)
    ann = AnnotationSet([model("G1", ex, g)])
    full = iso("full", "L1", ex, g, tpm=20)
    short5 = iso("short5", "L1", [(ex[0][0] + 100, ex[0][1])] + ex[1:], g, tpm=5)
    d, *_ = decide(ann, [full, short5], g, p=Params(ALT_END_ISOFORMS="all"))
    assert d["short5"].action == "new_isoform"                     # PASA: a UTR variant
    d, *_ = decide(ann, [full, short5], g)
    assert d["short5"].action == "alt_site_recorded" and d["short5"].status == 47


def test_an_end_without_utr_is_never_anchored(tmp_path):
    ex, _, g = _locus(tmp_path)
    m = model("G1", ex, g)
    t = m.transcripts["G1.t1"]
    cds_only = [(t.cds[0][0], ex[0][1])] + ex[1:]                   # 5' UTR trimmed away
    ann = AnnotationSet([model("G1", cds_only, g)])
    longer = iso("longer", "L1", ex, g)
    d, *_ = decide(ann, [longer], g, sites_for([longer]), Params(ADDITIVE=True))
    assert d["longer"].action == "utr_update"


def test_with_sites_a_model_whose_utr_holds_a_longer_orf_is_not_shortened(tmp_path):
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 300, 300])
    sense = [c for c in (a + b + d for a in "ACGT" for b in "ACGT" for d in "ACGT") if TABLE[c] not in "*M"]
    tail = "CCCATG" + "".join(L.rng.choice(sense) for _ in range(400)) + "TAACCC"
    L._add(tail)
    g = L.genome()
    m = model("G1", ex, g)
    t = m.transcripts["G1.t1"]
    t.exons = ex[:-1] + [(ex[-1][0], ex[-1][1] + len(tail))]
    ann = AnnotationSet([m])
    short = iso("short", "L1", ex, g)
    d, *_ = decide(ann, [short], g, p=Params(TRUST_FL_STATUS=True))     # no sites: PASA cuts it back
    assert d["short"].action == "utr_update"
    sites = sites_for([short], known=[("c", "+", TSS, ex[0][0])])
    d, decisions, splits, _ = decide(ann, [short], g, sites, Params(TRUST_FL_STATUS=True))
    assert d["short"].action == "matches_existing"
    new, _, _ = apply_round(ann, decisions, splits, "r1", {})
    assert new.genes["G1"].transcripts["G1.t1"].exons == t.exons
