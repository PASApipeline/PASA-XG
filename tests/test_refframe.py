"""Reference-frame check of isoform ORFs (annot/refframe.py)."""

from pasa_xg.annot import refframe
from pasa_xg.annot.checks import Params
from pasa_xg.isoforms import Orf
from pasa_xg.models import AnnotationSet

from pasa_xg.annot.structure import longest_orf

from test_annot import TABLE, Locus, iso, model


def _gene(tmp_path):
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 300, 300])
    g = L.genome()
    return ex, g, AnnotationSet([model("G1", ex, g)])


def test_isoform_starting_inside_the_cds_is_a_partial_in_the_reference_frame(tmp_path):
    ex, g, ann = _gene(tmp_path)
    ref = ann.genes["G1"].transcripts["G1.t1"]
    i = iso("frag", "L", ex[1:], g)                  # 5' end inside the CDS, no TSS support
    refframe.check({"frag": i}, ann, g, TABLE, None, Params())
    assert i.orf.orf_type == "5prime_partial" and not i.full_length
    assert i.tx.attrs["cds_start_NF"] == "true" and "cds_end_NF" not in i.tx.attrs
    assert g.protein(ref, TABLE).endswith(g.protein(i.tx, TABLE))   # the reference protein's tail


def test_a_supported_tss_keeps_the_isoforms_own_orf(tmp_path):
    ex, g, ann = _gene(tmp_path)
    i = iso("frag", "L", ex[1:], g)
    i.tss = True
    before = (i.orf, list(i.tx.cds))
    refframe.check({"frag": i}, ann, g, TABLE, None, Params())
    assert (i.orf, i.tx.cds) == before


def test_reference_start_opening_a_short_orf_makes_the_isoform_noncoding(tmp_path):
    ex, g, ann = _gene(tmp_path)
    shifted = [ex[0], (ex[1][0] + 1, ex[1][1]), ex[2]]   # one base lost: frame shift after exon 1
    i = iso("fs", "L", shifted, g)
    i.orf = Orf("fs.p1", "complete", 700, 790, 30)      # an ORF downstream of the reference ATG
    i.tx.cds = i.tx.segments_for(700, 790)
    refframe.check({"fs": i}, ann, g, TABLE, None, Params(), min_aa=500)
    assert i.orf is None and not i.tx.cds and i.tx.attrs["pasa_xg_orf"] == "reference_start_uorf"


def test_the_reference_structure_is_left_alone(tmp_path):
    ex, g, ann = _gene(tmp_path)
    i = iso("same", "L", ex, g)
    before = (i.orf, list(i.tx.cds))
    n = refframe.check({"same": i}, ann, g, TABLE, None, Params())
    assert (i.orf, i.tx.cds) == before and not n


def test_updated_models_are_checked_but_reference_cds_is_kept(tmp_path):
    ex, g, ann = _gene(tmp_path)
    new = ann.copy(label="u")
    ref = new.genes["G1"].transcripts["G1.t1"]
    frag = ref.copy(id="G1.pxg1", exons=ex[1:], attrs={})
    frag.cds = frag.segments_for(*longest_orf(g.transcript_seq(frag), TABLE, True, True))
    new.genes["G1"].transcripts["G1.pxg1"] = frag
    utr = (ex[2][0], ex[2][1] + 50)
    ref.exons = ex[:2] + [utr]                          # a UTR update: same CDS
    n = refframe.check_models(new, ann, g, TABLE, None, Params())
    assert ref.cds == ann.genes["G1"].transcripts["G1.t1"].cds
    assert frag.attrs.get("cds_start_NF") == "true" and n["partial in the reference frame"] == 1
