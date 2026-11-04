"""Noncoding pass (docs/NONCODING_DESIGN.md)."""

from pasa_xg.annot.checks import Params
from pasa_xg.annot.noncoding import run_noncoding
from pasa_xg.models import AnnotationSet

from test_annot import TABLE, Locus, iso, model


def _locus(tmp_path):
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 300, 300])
    L.spacer(2000)
    other = L.gene([300, 300, 300])          # not in the annotation: a novel locus
    g = L.genome()
    return ex, other, g


def _iso(name, exons, g, n, coding=True):
    i = iso(name, "L_" + name, exons, g)
    i.tx.attrs["n_samples"] = str(n)
    if not coding:
        i.tx.cds = []
    return i


def run(ann, isos, g):
    new, d = run_noncoding(ann, {i.id: i for i in isos}, [i.id for i in isos], g, TABLE, Params(), None)
    return new, {x.isoform_id: x for x in d}


def test_retained_intron_isoform_is_added_as_noncoding_with_enough_samples(tmp_path):
    ex, _, g = _locus(tmp_path)
    ann = AnnotationSet([model("G1", ex, g)])
    retained = [ex[0], (ex[1][0], ex[2][1])]
    new, d = run(ann, [_iso("ri", retained, g, 3, coding=False)], g)
    assert d["ri"].action == "noncoding_isoform" and d["ri"].annotation_id == "G1.pxg1"
    t = new.genes["G1"].transcripts["G1.pxg1"]
    assert not t.cds and t.attrs["pasa_xg_coding"] == "none" and t.exons == retained
    assert new.genes["G1"].attrs["pasa_xg_gene_class"] == "coding"
    _, d = run(ann, [_iso("ri", retained, g, 2, coding=False)], g)      # 2 samples, no supported end
    assert d["ri"].status == 103


def test_isoform_that_assembles_with_a_transcript_adds_nothing(tmp_path):
    ex, _, g = _locus(tmp_path)
    ann = AnnotationSet([model("G1", ex, g)])
    new, d = run(ann, [_iso("frag", ex[1:], g, 5, coding=False)], g)
    assert d["frag"].status == 104 and len(new.genes["G1"].transcripts) == 1


def test_utr_variant_with_the_gene_protein_is_low_coding_at_the_lower_bar(tmp_path):
    ex, _, g = _locus(tmp_path)
    ann = AnnotationSet([model("G1", ex, g)])
    utr_intron = ex[:2] + [(ex[2][0], ex[2][0] + 250), (ex[2][0] + 280, ex[2][1])]
    new, d = run(ann, [_iso("utr", utr_intron, g, 2)], g)
    assert d["utr"].action == "noncoding_isoform"
    t = new.genes["G1"].transcripts[d["utr"].annotation_id]
    assert t.attrs["pasa_xg_coding"] == "low" and t.cds
    _, d = run(ann, [_iso("utr", utr_intron, g, 1)], g)                  # 1 sample, no supported end
    assert d["utr"].status == 103


def test_spliced_isoform_at_an_unannotated_locus_is_a_novel_noncoding_gene(tmp_path):
    ex, other, g = _locus(tmp_path)
    ann = AnnotationSet([model("G1", ex, g)])
    new, d = run(ann, [_iso("nov", other, g, 3, coding=False), _iso("nov2", other[1:], g, 3, coding=False)], g)
    assert d["nov"].action == "noncoding_gene" and d["nov2"].status == 104
    gene = new.genes[d["nov"].genes[0]]
    assert gene.attrs["pasa_xg_gene_class"] == "likely_noncoding" and gene.attrs["pasa_xg_novel"] == "noncoding"
    _, d = run(ann, [_iso("mono", [other[0]], g, 9, coding=False)], g)  # mono-exonic, no gene: left out
    assert d["mono"].status == 106


def test_an_embedded_small_rna_gene_does_not_make_an_isoform_multigenic(tmp_path):
    from pasa_xg.models import Gene, Transcript
    ex, _, g = _locus(tmp_path)
    mir = Transcript("MIR.t1", "MIR", "c", "+", [(ex[1][1] + 50, ex[1][1] + 130)])     # in G1's 2nd intron
    ann = AnnotationSet([model("G1", ex, g), Gene("MIR", "c", "+", {"MIR.t1": mir})])
    retained = [ex[0], (ex[1][0], ex[2][1])]
    _, d = run(ann, [_iso("ri", retained, g, 3, coding=False)], g)
    assert d["ri"].action == "noncoding_isoform" and d["ri"].genes == ["G1"]


def test_isoform_of_a_gene_with_a_readthrough_goes_to_the_better_matching_gene(tmp_path):
    from pasa_xg.models import Gene, Transcript
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 300, 300, 300])
    L.spacer(500)
    ex2 = L.gene([300, 300])
    g = L.genome()
    rt = Transcript("RT.t1", "RT", "c", "+", ex[:-1] + [(ex[-1][0], ex[-1][1])] + ex2)   # readthrough gene
    ann = AnnotationSet([model("G1", ex, g), Gene("RT", "c", "+", {"RT.t1": rt})])
    skip = [ex[0], ex[2], ex[3]]                                           # exon skip within G1
    _, d = run(ann, [_iso("sk", skip, g, 3, coding=False)], g)
    assert d["sk"].genes == ["G1"]


def test_without_sample_counts_every_isoform_is_evaluated_as_it_is(tmp_path):
    ex, other, g = _locus(tmp_path)
    ann = AnnotationSet([model("G1", ex, g)])
    retained = _iso("ri", [ex[0], (ex[1][0], ex[2][1])], g, 1, coding=False)
    novel = _iso("nov", other, g, 1, coding=False)
    for i in (retained, novel):
        del i.tx.attrs["n_samples"]                                     # a single-sample run
    _, d = run(ann, [retained, novel], g)
    assert d["ri"].action == "noncoding_isoform" and d["nov"].action == "noncoding_gene"


def test_isoform_matching_two_genes_on_shared_introns_is_flagged_and_skipped(tmp_path):
    from pasa_xg.models import Gene, Transcript
    ex, _, g = _locus(tmp_path)
    dup = Transcript("G2.t1", "G2", "c", "+",                            # a second gene on G1's 2nd intron
                     [ex[1], (ex[2][0], ex[2][0] + 100), (ex[2][0] + 150, ex[2][1])])
    ann = AnnotationSet([model("G1", ex, g), Gene("G2", "c", "+", {"G2.t1": dup})])
    retained = [(ex[0][0], ex[1][1]), ex[2]]
    new, d = run(ann, [_iso("ri", retained, g, 3, coding=False)], g)
    assert d["ri"].status == 105 and sorted(d["ri"].genes) == ["G1", "G2"]
    assert d["ri"].comment.startswith("shares introns") and len(new.genes["G1"].transcripts) == 1
