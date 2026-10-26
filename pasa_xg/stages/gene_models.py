"""Stage 6: coding gene models straight from the isoforms."""

from __future__ import annotations

import os
from typing import Any, Dict

from .. import isoforms as iso_mod
from ..config import Settings
from ..genome import Genome, codon_table, translate
from ..gff import write_bed, write_fasta, write_gff3
from ..models import AnnotationSet, Gene
from ..runner import Run, code_digest


def build(isos: Dict[str, iso_mod.Isoform]) -> AnnotationSet:
    """Genes from LRAA's gene grouping; one transcript per isoform with an ORF."""
    genes: Dict[str, Gene] = {}
    for iso in isos.values():
        if not iso.orf:
            continue
        t = iso.tx.copy()
        t.attrs = {"orf_type": iso.orf.orf_type, "orf_aa_len": str(iso.orf.aa_len),
                   "TPM": f"{iso.tpm:.3f}", "isoform_fraction": f"{iso.isoform_fraction:.3f}"}
        g = genes.get(t.gene_id)
        if g is None:
            g = genes[t.gene_id] = Gene(t.gene_id, t.contig, t.strand)
        g.add(t)
    return AnnotationSet(genes.values(), label="gene_models")


def run(run: Run, s: Settings, prep: Dict[str, Any], coding: Dict[str, Any]) -> Dict[str, Any]:
    inputs = {"isoforms": coding["isoforms_path"], "code": s.genetic_code, "stage": "gene_models/1", "src": code_digest("pasa_xg.stages.gene_models", "pasa_xg.gff", "pasa_xg.models")}

    def body(d: str) -> Dict[str, Any]:
        isos = iso_mod.load(coding["isoforms_path"])
        models = build(isos)
        genome = Genome(prep["genome_path"])
        table = codon_table(s.genetic_code)
        gff = os.path.join(d, "gene_models.gff3")
        with open(gff, "w") as fh:
            write_gff3(models, fh)
        with open(os.path.join(d, "gene_models.bed"), "w") as fh:
            write_bed(models.transcripts(), fh)
        with open(os.path.join(d, "gene_models.pep"), "w") as fh:
            write_fasta(((t.id, translate(genome.cds_seq(t), table)) for t in models.transcripts()), fh)
        with open(os.path.join(d, "gene_models.cds"), "w") as fh:
            write_fasta(((t.id, genome.cds_seq(t)) for t in models.transcripts()), fh)
        n_tx = sum(len(g.transcripts) for g in models)
        return {"gff3_path": gff, "genes": len(models), "transcripts": n_tx}

    return run.run_stage("gene_models", "06_gene_models", inputs, body)
