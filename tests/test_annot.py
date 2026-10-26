"""Annotation update: structure operations, and one synthetic locus per decision."""

import random

import pytest

from pasa_xg.annot.apply import apply_round
from pasa_xg.annot.checks import Params, align_proteins, utr_exon_count
from pasa_xg.annot.compare import run_round
from pasa_xg.annot.structure import longest_orf, merge_compatible, stitch
from pasa_xg.genome import Genome, codon_table
from pasa_xg.isoforms import Isoform, Orf
from pasa_xg.models import AnnotationSet, Gene, Transcript

TABLE = codon_table()


# ------------------------------------------------------------------ structure ops

def test_merge_compatible_extends_and_rejects_conflicts():
    a = [(100, 200), (300, 400), (500, 600)]
    b = [(350, 400), (500, 700)]
    assert merge_compatible(a, b) == [(100, 200), (300, 400), (500, 700)]
    assert merge_compatible(a, [(150, 450)]) is None                 # covers an intron
    assert merge_compatible(a, [(350, 400), (520, 600)]) is None     # different acceptor
    assert merge_compatible(a, [(390, 410)]) == a                     # 10 bp into an intron: fuzz
    assert merge_compatible(a, [(390, 430)]) is None                  # 30 bp: beyond fuzz
    assert merge_compatible(a, [(800, 900)]) is None                  # no overlap


def test_stitch_replaces_structure_across_the_thread():
    template = [(100, 200), (300, 400), (500, 600), (700, 800)]
    thread = [(350, 400), (450, 480), (500, 560)]   # adds an exon between 400 and 500
    assert stitch(template, thread) == [(100, 200), (300, 400), (450, 480), (500, 600), (700, 800)]
    # a single-exon thread spanning an intron makes it a retained intron
    assert stitch(template, [(380, 520)]) == [(100, 200), (300, 600), (700, 800)]


def test_longest_orf_complete_and_partial():
    seq = "CC" + "ATG" + "GCT" * 10 + "TAA" + "CC"
    assert longest_orf(seq, TABLE) == (3, 38)
    open_end = "C" + "GCT" * 5 + "ATG" + "GCT" * 4
    assert longest_orf(open_end, TABLE) is None
    assert longest_orf(open_end, TABLE, allow_3partial=True) is not None
    assert longest_orf("GCT" * 5 + "TAA", TABLE, allow_5partial=True) == (1, 18)


def test_protein_alignment_and_utr_exons():
    pid, length = align_proteins("MKVLAAGIVLLAQW*", "MKVLAAGIVLLAQW*")
    assert pid == 100 and length == 14
    t = Transcript("t", "g", "c", "+", [(1, 10), (20, 30), (40, 50), (60, 70)], cds=[(60, 65)])
    assert utr_exon_count(t) == 3


# ------------------------------------------------------------------ synthetic loci

class Locus:
    """A genome built from gene specs: each gene is a list of exon lengths and an ORF."""

    def __init__(self, tmp_path, seed=3):
        self.rng = random.Random(seed)
        self.seq = []
        self.pos = 0
        self.tmp = tmp_path

    def _add(self, s):
        self.seq.append(s)
        self.pos += len(s)

    def spacer(self, n):
        self._add("".join(self.rng.choice("CT") for _ in range(n)))

    def gene(self, exon_lens, utr5=60, intron=200):
        """Lay down a gene; returns its exons. The transcript is UTR (C/T only, so no
        ATG or stop) + ATG + non-stop codons + TAA + UTR."""
        total = sum(exon_lens)
        n_codons = (total - utr5 - 60 - 6) // 3
        sense = [c for c in (a + b + d for a in "ACGT" for b in "ACGT" for d in "ACGT")
                 if TABLE[c] not in "*M"]
        cds = "ATG" + "".join(self.rng.choice(sense) for _ in range(n_codons)) + "TAA"
        tx = "".join(self.rng.choice("CT") for _ in range(utr5)) + cds
        tx += "".join(self.rng.choice("CT") for _ in range(total - len(tx)))
        exons, k = [], 0
        for i, n in enumerate(exon_lens):
            exons.append((self.pos + 1, self.pos + n))
            self._add(tx[k:k + n])
            k += n
            if i < len(exon_lens) - 1:
                self._add("GT" + "".join(self.rng.choice("CT") for _ in range(intron - 4)) + "AG")
        return exons

    def genome(self):
        self.spacer(500)
        p = self.tmp / "g.fa"
        p.write_text(">c\n" + "".join(self.seq) + "\n")
        return Genome(str(p))


def model(gid, exons, genome, strand="+"):
    t = Transcript(f"{gid}.t1", gid, "c", strand, exons)
    orf = longest_orf(genome.transcript_seq(t), TABLE)
    t.cds = t.segments_for(*orf)
    return Gene(gid, "c", strand, {t.id: t})


def iso(iid, lraa_gene, exons, genome, strand="+", tpm=10.0):
    """An isoform with its best ORF: complete if there is one, else partial."""
    t = Transcript(iid, lraa_gene, "c", strand, exons)
    seq = genome.transcript_seq(t)
    orf, kind = longest_orf(seq, TABLE), "complete"
    if orf is None:
        orf, kind = longest_orf(seq, TABLE, True, True), "internal"
    o = None
    if orf:
        o = Orf(iid + ".p1", kind, orf[0], orf[1], (orf[1] - orf[0] + 1) // 3)
        t.cds = t.segments_for(*orf)
    return Isoform(t, tpm=tpm, isoform_fraction=1.0, orf=o)


def decide(ann, isos, genome, p=None):
    decisions, splits, _ = run_round(ann, isos, genome, TABLE, p or Params())
    return {d.isoform_id: d for d in decisions}, decisions, splits


def test_identical_isoform_matches_and_longer_utr_updates(tmp_path):
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 300, 300])
    g = L.genome()
    ann = AnnotationSet([model("G1", ex, g)])
    same = iso("same", "L1", ex, g)
    d, _, _ = decide(ann, [same], g)
    assert d["same"].action == "matches_existing"
    shorter_model = AnnotationSet([model("G1", [(ex[0][0] + 30, ex[0][1])] + ex[1:], g)])
    d, decisions, splits = decide(shorter_model, [same], g)
    assert d["same"].action == "utr_update"
    new, applied, changed = apply_round(shorter_model, decisions, splits, "r1", {})
    assert new.genes["G1"].transcripts["G1.t1"].exons == ex and changed == {"G1"}


def test_exon_skip_isoform_becomes_alternative_isoform(tmp_path):
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 150, 150, 300])
    g = L.genome()
    ann = AnnotationSet([model("G1", ex, g)])
    full = iso("full", "L1", ex, g)
    skip = iso("skip", "L1", [ex[0], ex[2], ex[3]], g)  # 150 bp: frame kept
    d, decisions, splits = decide(ann, [full, skip], g)
    assert d["full"].action == "matches_existing"
    assert d["skip"].action == "new_isoform" and d["skip"].status == 9
    new, applied, _ = apply_round(ann, decisions, splits, "r1", {})
    assert set(new.genes["G1"].transcripts) == {"G1.t1", "G1.pxg1"}
    assert applied["skip"] == "G1.pxg1"


def test_intergenic_isoform_is_a_novel_gene(tmp_path):
    L = Locus(tmp_path)
    L.spacer(300)
    ex1 = L.gene([300, 300, 300])
    L.spacer(3000)
    ex2 = L.gene([300, 300, 300])
    g = L.genome()
    ann = AnnotationSet([model("G1", ex1, g)])
    d, decisions, splits = decide(ann, [iso("nov", "L2", ex2, g)], g)
    assert d["nov"].action == "novel_gene"
    new, applied, _ = apply_round(ann, decisions, splits, "r1", {})
    assert applied["nov"] == "PXG_gene_1.t1" and "PXG_gene_1" in new.genes


def test_isoform_spanning_two_genes_merges_them(tmp_path):
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 300, 300, 300])
    g = L.genome()
    # the annotation splits one real gene in two
    ann = AnnotationSet([model("GA", ex[:2], g), model("GB", ex[2:], g)])
    # GA, GB lack proper ORFs as fragments; give them their own CDS from the full ORF
    full_orf_t = Transcript("x", "x", "c", "+", ex)
    orf = longest_orf(g.transcript_seq(full_orf_t), TABLE)
    cds = full_orf_t.segments_for(*orf)
    for gid in ("GA", "GB"):
        t = next(iter(ann.genes[gid].transcripts.values()))
        t.cds = [(max(a, t.start), min(b, t.end)) for a, b in cds if a <= t.end and b >= t.start]
    d, decisions, splits = decide(ann, [iso("span", "L1", ex, g)], g)
    assert d["span"].action == "merge_genes"
    new, applied, _ = apply_round(ann, decisions, splits, "r1", {})
    assert set(new.genes) == {"GA"} and new.genes["GA"].attrs["merged_from"] == "GA,GB"


def test_opposite_strand_isoform_is_antisense(tmp_path):
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 300, 300])
    g = L.genome()
    ann = AnnotationSet([model("G1", ex, g)])
    # an unspliced, noncoding transcript on the other strand inside the gene
    anti = Isoform(Transcript("anti", "L9", "c", "-", [(ex[1][0], ex[1][1])]), tpm=5)
    d, _, _ = decide(ann, [anti], g)
    assert d["anti"].action == "antisense"


def test_second_round_places_an_isoform_once_its_genes_are_merged(tmp_path):
    from pasa_xg.config import Settings
    from pasa_xg.stages.annot_update import update
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 150, 150, 300, 300])
    g = L.genome()
    full_t = Transcript("x", "x", "c", "+", ex)
    cds = full_t.segments_for(*longest_orf(g.transcript_seq(full_t), TABLE))
    genes = []
    for gid, part in (("GA", ex[:2]), ("GB", ex[2:])):
        t = Transcript(f"{gid}.t1", gid, "c", "+", part,
                       [(max(a, part[0][0]), min(b, part[-1][1])) for a, b in cds
                        if a <= part[-1][1] and b >= part[0][0]])
        genes.append(Gene(gid, "c", "+", {t.id: t}))
    ann = AnnotationSet(genes)
    isos = {"span": iso("span", "L1", ex, g),
            "skip": iso("skip", "L1", [ex[0], ex[2], ex[3], ex[4]], g)}
    s = Settings(genome=str(tmp_path / "g.fa"), input_type="transcripts", out_dir=str(tmp_path),
                 reads=["r.fa"])
    out = update(ann, isos, g, s, str(tmp_path / "upd"), Params())
    import pandas as pd
    d = pd.read_csv(tmp_path / "upd" / "annot_update.all_rounds.tsv", sep="\t")
    r1 = d[d["round"] == 1].set_index("isoform_id")
    r2 = d[d["round"] == 2].set_index("isoform_id")
    assert r1.loc["span", "action"] == "merge_genes"
    assert r1.loc["skip", "action"] == "rejected"
    assert r2.loc["skip", "action"] == "new_isoform"
    assert out["incorporated"] == 2


def test_nmd_distance_stop_to_last_junction():
    from pasa_xg.annot.checks import nmd_distance
    # exons 100 + 100 + 100 bases; stop codon ends at transcript base 30 (first exon)
    t = Transcript("t", "g", "c", "+", [(1, 100), (201, 300), (401, 500)], cds=[(10, 30)])
    assert nmd_distance(t) == 200 - 30
    t2 = t.copy(cds=[(10, 100), (201, 300), (401, 420)])  # stop in the last exon
    assert nmd_distance(t2) < 0
    m = Transcript("m", "g", "c", "-", [(1, 100), (201, 300), (401, 500)], cds=[(470, 490)])
    assert nmd_distance(m) == 200 - 31   # minus strand: stop ends 31 bases into the transcript
    assert nmd_distance(Transcript("s", "g", "c", "+", [(1, 500)], cds=[(10, 30)])) == 0


def test_update_history_accumulates_across_chained_runs(tmp_path):
    import io as _io
    from pasa_xg.gff import read_annotation, write_gff3
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 300, 300])
    g = L.genome()
    trimmed = [(ex[0][0] + 60, ex[0][1])] + ex[1:-1] + [(ex[-1][0], ex[-1][1] - 60)]
    start = AnnotationSet([model("G1", trimmed, g)])
    # run 1 restores the 5' end only
    five = iso("five", "L1", ex[:-1] + [(ex[-1][0], ex[-1][1] - 60)], g)
    _, decisions, splits = decide(start, [five], g)
    run1, _, _ = apply_round(start, decisions, splits, "round_01", {}, run_label="gen2")
    # hand the output to the next run through GFF3, as a chained update would
    out = _io.StringIO()
    write_gff3(run1, out)
    p = tmp_path / "run1.gff3"
    p.write_text(out.getvalue())
    reloaded = read_annotation(str(p))
    full = iso("full", "L1", ex, g)
    _, decisions, splits = decide(reloaded, [full], g)
    run2, _, _ = apply_round(reloaded, decisions, splits, "round_01", {}, run_label="gen3")
    t = run2.genes["G1"].transcripts["G1.t1"]
    assert t.exons == ex
    assert t.attrs["pasa_xg_history"] == "gen2/round_01:utr_update:five,gen3/round_01:utr_update:full"
    assert t.attrs["pasa_xg_run"] == "gen3"


def test_additive_mode_keeps_an_unsupported_model(tmp_path):
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 150, 150, 300])
    g = L.genome()
    skip_model = AnnotationSet([model("G1", [ex[0], ex[2], ex[3]], g)])  # no read shows it
    full = iso("full", "L1", ex, g)
    d, decisions, splits = decide(skip_model, [full], g)
    assert d["full"].action == "protein_update" and d["full"].status == 8   # PASA: replaced
    d, decisions, splits = decide(skip_model, [full], g, Params(ADDITIVE=True))
    assert d["full"].action == "new_isoform"
    new, _, _ = apply_round(skip_model, decisions, splits, "r1", {})
    tx = new.genes["G1"].transcripts
    assert tx["G1.t1"].exons == [ex[0], ex[2], ex[3]] and len(tx) == 2
