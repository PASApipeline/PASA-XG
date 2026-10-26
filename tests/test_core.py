import io
import json
import os

import pytest

from pasa_xg.genome import Genome, codon_table, revcomp, translate
from pasa_xg.gff import read_annotation, write_gff3, write_gtf
from pasa_xg.models import AnnotationSet, Gene, Transcript
from pasa_xg.runner import Run
from pasa_xg.stages.coding import parse_pep_headers


def test_transcript_geometry_plus_and_minus():
    t = Transcript("t", "g", "c", "+", [(100, 110), (200, 210)], cds=[(105, 110), (200, 205)])
    assert t.introns() == [(111, 199)]
    assert t.length() == 22
    assert t.transcript_to_genome(1) == 100 and t.transcript_to_genome(12) == 200
    assert t.segments_for(6, 17) == [(105, 110), (200, 205)]
    assert t.cds_phases() == [0, 0]
    assert t.utrs() == ([(100, 104)], [(206, 210)])
    m = Transcript("m", "g", "c", "-", [(100, 110), (200, 210)])
    assert m.transcript_to_genome(1) == 210 and m.transcript_to_genome(12) == 110
    assert m.segments_for(1, 13) == [(109, 110), (200, 210)]


def test_cds_phases_track_split_codons():
    t = Transcript("t", "g", "c", "+", [(1, 10), (21, 40)], cds=[(3, 10), (21, 30)])
    assert t.cds_phases() == [0, 1]  # 8 bases in the first segment leave 1 to finish a codon
    m = Transcript("m", "g", "c", "-", [(1, 10), (21, 40)], cds=[(3, 10), (21, 30)])
    assert m.cds_phases() == [2, 0]  # minus strand: the 10 bases at 21..30 come first


def test_translate_and_revcomp():
    tab = codon_table()
    assert translate("ATGTTTTAA", tab) == "MF*"
    assert revcomp("ATGC") == "GCAT"
    assert translate("TGA", codon_table("Mitochondrial-Vertebrates")) == "W"


GTF = """c1\tsrc\ttranscript\t11\t60\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
c1\tsrc\texon\t11\t20\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
c1\tsrc\texon\t41\t60\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
c1\tsrc\tCDS\t14\t20\t.\t+\t0\tgene_id "G1"; transcript_id "T1";
c1\tsrc\tCDS\t41\t49\t.\t+\t2\tgene_id "G1"; transcript_id "T1";
c1\tsrc\tstop_codon\t50\t52\t.\t+\t0\tgene_id "G1"; transcript_id "T1";
"""


def test_gtf_stop_codon_joins_cds_and_round_trips(tmp_path):
    p = tmp_path / "a.gtf"
    p.write_text(GTF)
    ann = read_annotation(str(p))
    t = next(ann.transcripts())
    assert t.exons == [(11, 20), (41, 60)] and t.cds == [(14, 20), (41, 52)]
    out = io.StringIO()
    write_gff3(ann, out)
    p3 = tmp_path / "a.gff3"
    p3.write_text(out.getvalue())
    t3 = next(read_annotation(str(p3)).transcripts())
    assert (t3.id, t3.gene_id, t3.exons, t3.cds) == ("T1", "G1", t.exons, t.cds)
    out = io.StringIO()
    write_gtf(ann, out)
    p2 = tmp_path / "b.gtf"
    p2.write_text(out.getvalue())
    assert next(read_annotation(str(p2)).transcripts()).cds == t.cds


def test_overlap_index():
    genes = [Gene(f"g{i}", "c", "+") for i in range(3)]
    for i, g in enumerate(genes):
        g.add(Transcript(f"t{i}", g.id, "c", "+", [(i * 1000 + 1, i * 1000 + 500)]))
    big = Gene("big", "c", "-")
    big.add(Transcript("tb", "big", "c", "-", [(1, 5000)]))
    ann = AnnotationSet(genes + [big])
    assert {g.id for g in ann.overlapping("c", 1400, 2100)} == {"g1", "g2", "big"}
    assert {g.id for g in ann.overlapping("c", 1400, 2100, "+")} == {"g1", "g2"}
    assert ann.overlapping("c", 600, 900, "+") == []


def test_genome_spliced_minus(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">c\nAAAACCCCGGGGTTTT\n")
    g = Genome(str(fa))
    assert g.spliced("c", "+", [(1, 2), (9, 10)]) == "AAGG"
    assert g.spliced("c", "-", [(1, 2), (9, 10)]) == "CCTT"


def test_parse_transdecoder_headers(tmp_path):
    pep = tmp_path / "x.pep"
    pep.write_text(
        ">t:gi|1:+:comp-1:iso-1.p1 GENE.x~~t.p1  ORF type:complete (+),score=40.1 len:120 "
        "t:gi|1:+:comp-1:iso-1:11-373(+)\nMAAA\n"
        ">t2.p1 GENE.y~~t2.p1  ORF type:5prime_partial (-),score=3 len:101 t2:1-303(-)\nAAA\n")
    orfs = parse_pep_headers(str(pep))
    assert list(orfs) == ["t:gi|1:+:comp-1:iso-1"]
    o = orfs["t:gi|1:+:comp-1:iso-1"]
    assert (o.orf_type, o.t_start, o.t_end, o.aa_len, o.complete) == ("complete", 11, 373, 120, True)


def test_stage_checkpoint_skip_rerun_and_force(tmp_path):
    src = tmp_path / "in.txt"
    src.write_text("one")
    calls = []

    def body(d):
        calls.append(1)
        p = os.path.join(d, "o.txt")
        open(p, "w").write(src.read_text())
        return {"o_path": p, "n": 1}

    def go(force=None):
        r = Run(str(tmp_path / "run"), {}, force_from=force)
        out = r.run_stage("s1", "s1", {"in": str(src)}, body)
        r.run_stage("s2", "s2", {"up": out["o_path"]}, body)
        return out

    go()
    assert len(calls) == 2
    out = go()
    assert len(calls) == 2 and out["o_path"].endswith("s1/o.txt") and out["n"] == 1
    src.write_text("two")  # input content changed: both stages rerun
    go()
    assert len(calls) == 4
    go(force="s1")  # forcing a stage reruns everything after it
    assert len(calls) == 6
    m = json.load(open(tmp_path / "run" / "manifest.json"))
    assert m["stages"]["s2"]["status"] == "done"


def test_five_prime_incomplete_cds_phase_is_honored(tmp_path):
    # CDS starts 1 base into a codon (phase 1): the protein must start at the next base
    fa = tmp_path / "g.fa"
    fa.write_text(">c\n" + "C" * 10 + "A" + "ATGTTTTAA" + "C" * 10 + "\n")
    gtf = tmp_path / "a.gtf"
    gtf.write_text('c\ts\texon\t1\t30\t.\t+\t.\tgene_id "G"; transcript_id "T";\n'
                   'c\ts\tCDS\t11\t20\t.\t+\t1\tgene_id "G"; transcript_id "T";\n')
    t = next(read_annotation(str(gtf)).transcripts())
    assert t.cds_offset == 1 and t.cds_phases() == [1]
    assert Genome(str(fa)).protein(t, codon_table()) == "MF*"
    out = io.StringIO()
    write_gff3(read_annotation(str(gtf)), out)
    assert "\tCDS\t11\t20\t.\t+\t1\t" in out.getvalue()
