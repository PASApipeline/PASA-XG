"""GENCODE attributes through read and write, and the guards on protein updates."""

import io

from pasa_xg.annot.checks import Params, compare_proteins
from pasa_xg.gff import read_annotation, write_gtf

from test_annot import TABLE, Locus, model


def test_repeated_tags_survive_read_and_write(tmp_path):
    p = tmp_path / "a.gtf"
    p.write_text('c\tX\ttranscript\t1\t100\t.\t+\t.\tgene_id "G"; transcript_id "T"; gene_name "N"; tag "basic"; tag "MANE_Select";\n'
                 'c\tX\texon\t1\t100\t.\t+\t.\tgene_id "G"; transcript_id "T";\n')
    t = next(iter(read_annotation(str(p)).transcripts()))
    assert t.attrs["tag"] == "basic,MANE_Select"
    out = io.StringIO()
    write_gtf(read_annotation(str(p)), out)
    line = out.getvalue().splitlines()[0]
    assert 'tag "basic";' in line and 'tag "MANE_Select";' in line and 'gene_name "N";' in line


def test_frozen_and_complete_cds_guards(tmp_path):
    L = Locus(tmp_path)
    L.spacer(300)
    ex = L.gene([300, 300, 300])
    g = L.genome()
    old = model("G1", ex, g).transcripts["G1.t1"]
    shorter = old.copy(cds=[old.cds[0], old.cds[1]] if len(old.cds) > 2 else old.cds)
    shorter.cds = shorter.segments_for(61, 600)              # same start, no stop codon: a partial protein
    ok, why = compare_proteins(old, shorter, g, TABLE, Params())
    assert not ok and "partial" in why
    old.attrs = {"tag": "basic,seleno"}
    ext = old.copy(exons=old.exons, cds=old.cds)
    ext.cds = ext.segments_for(61, 663)
    ok, why = compare_proteins(old, ext, g, TABLE, Params())
    assert not ok and "seleno" in why
