# PASA-XG design

Draft, 2026-10-25. For review before the repo (`PASApipeline/PASA-XG`) is created.

PASA-XG ("cross-generation") is a rebuild of the PASA pipeline for every generation
of transcript evidence: ESTs and cDNAs, short-read assemblies (Trinity, StringTie)
and PacBio/ONT long reads. It keeps what PASA is used for (updating a genome
annotation from transcript evidence, gene models from transcripts, alternative
splicing) and replaces the engine: minimap2 aligns, LRAA assembles and quantifies,
TransDecoder finds the coding regions.

It does not reproduce Perl PASA's results. PyPASA
(`PASApipeline/PyPASA`, branch `perl-parity-in-memory`) remains the record of the
exact port of PASA's alignment assembly.

## 1. Scope

**v1:**

1. Align reads or transcripts with minimap2, or start from a BAM.
2. Assemble and quantify isoforms with LRAA (ref-guided when an annotation is given).
3. Find coding regions with TransDecoder; map them to the genome.
4. Compare the isoforms to the annotation and update it, over several rounds.
5. Gene models from the isoforms, with or without an annotation.
6. Alternative splicing classification.

**Not planned:** PASA's comprehensive transcriptome (merging de novo and
genome-guided assemblies), the MySQL/SQLite database, and the Perl web portal. A
new results browser may come later, built on the run outputs.

**v1 limits:** one input type per run (reads or transcripts, not both); minimap2 is
the only aligner.

## 2. Inputs

| Input | Required | Notes |
|---|---|---|
| Genome FASTA | yes | |
| Reads or transcripts (FASTA/FASTQ, gz ok), **or** a BAM | one of them | A BAM skips alignment. It must be coordinate-sorted and indexed (indexed by PASA-XG if not). |
| `--input_type` | yes | See below. |
| Annotation (GTF or GFF3) | no | Used for ref-guided LRAA and for the annotation update. The update needs CDS features; an annotation without CDS still guides LRAA, and the update is skipped with a message. |
| Pfam HMMs / protein DB | no | Homology evidence for TransDecoder. |

`--input_type` sets the aligner and LRAA options together; each can still be set
individually.

| `--input_type` | minimap2 | LRAA |
|---|---|---|
| `pacbio_hifi` | `-ax splice:hq -ub` | `--HiFi` |
| `ont` | `-ax splice -ub` | defaults |
| `transcripts` (Trinity, StringTie, ESTs, cDNAs) | PASA's current settings: `-ax splice --secondary=no -O6,24 -B4 -L -cs -ub` | `--pasa_mode --preassembled_input` |

`--minimap2_opts` overrides the minimap2 settings. For reads, `--chain_fragments`
adds LRAA's `--pasa_mode` (chaining reads that share an intron), for fragmented or
degraded samples.

**Strand.** No strand option: minimap2 always runs with `-ub`, so it picks each
spliced read's transcribed strand from its splice motifs and records it in the
`ts` tag, and LRAA resolves strands from `ts` on its own (a read moves to the
other strand only when every junction has a canonical motif on that strand).
Stranded and unstranded libraries are handled the same way. Unspliced reads carry
no strand evidence and stay on the strand they aligned to, which is right for
stranded libraries. A BAM given as input should come from minimap2 (or carry `ts`
tags). `--minimap2_opts` can still force `-uf`.

## 3. Pipeline

```
genome ─┬─ prep ── align (minimap2) ── assemble (LRAA) ── isoforms ── coding (TransDecoder)
annot ──┘                                  ▲                               │
          (ref-guided when annotation) ────┘           ┌───────────────────┼──────────────────┐
                                                annot update (rounds)  gene models     alt splicing
                                                       └───────────────────┴──────────────────┘
                                                                        report
```

| # | Stage | Does | Main outputs |
|---|---|---|---|
| 0 | `prep` | Index the genome (faidx; minimap2 `.mmi`, cacheable across runs). Load and check the annotation: parse GTF/GFF3 into gene models, verify CDS (translates, start/stop, frame), and report problems. | `genome.fa.fai`, `genome.mmi`, `annotation.models.pkl`, `annotation_check.tsv` |
| 1 | `align` | minimap2 → sorted, indexed BAM. Skipped for BAM input. | `aligned.bam` |
| 2 | `assemble` | LRAA, ref-guided if an annotation is given. | LRAA GTF, `quant.expr`, `quant.tracking`, TSS/PolyA beds |
| 3 | `isoforms` | Load the LRAA outputs into one isoform table: structure, gene, expression (reads, TPM, isoform fraction), TSS/polyA support, ref match class. Write isoform sequences. | `isoforms.parquet`, `isoforms.fa` |
| 4 | `coding` | TransDecoder LongOrfs + Predict on the isoform sequences (optional hmmscan / BLASTP or DIAMOND evidence); map the ORFs to genome coordinates. | `isoform_orfs.parquet`, `isoforms.coding.gff3` |
| 5 | `annot_update` | Compare isoforms to the annotation and update it, in rounds (§5). Needs an annotation with CDS. | `rounds/round_NN/`, `annotation.updated.gff3`, `annot_update.decisions.tsv` |
| 6 | `gene_models` | Gene models straight from the isoforms (§6). | `gene_models.gff3` |
| 7 | `alt_splice` | Alternative splicing events between isoforms (§7). | `alt_splice.events.tsv` |
| 8 | `report` | Summary tables, final outputs collected in one place. | `outputs/`, `summary.tsv` |

Stages 5 to 7 are independent of each other; each reads stages 3 and 4.

## 4. Data model

Plain Python dataclasses, held in memory during a stage and written as checkpoints
between stages. Coordinates are 1-based, inclusive, as in GFF.

- `Exon(start, end)`, `CDSSegment(start, end, phase)`.
- `Transcript`: id, gene_id, contig, strand, exons, CDS segments (optional),
  attributes. One class serves reference transcripts, LRAA isoforms and updated
  models.
- `Gene`: id, transcripts, attributes.
- `Isoform` (LRAA): a `Transcript` plus expression (`uniq_reads`, `all_reads`, TPM,
  `isoform_fraction`), TSS/polyA support, LRAA's ref match class, and its ORF from
  stage 4 (type: complete, 5'partial, 3'partial, internal; or none).
- `AnnotationSet`: genes by id with an interval index per contig/strand; a version
  label (`input`, `round_01`, ...).
- `Decision`: isoform id, round, action (§5.2), target gene(s), validation results,
  reason when rejected.

The `Transcript` / `Gene` model replaces the parts of Perl's `Gene_obj.pm`
(5,500 lines) that the update logic needs: CDS from exons and ORF, spliced
sequence, translation with a selectable genetic code, UTR derivation, structure
comparison.

## 5. Annotation update

### 5.1 One round

A port of PASA's `cDNA_annotation_comparer.dbi` (3,500 lines) and its helpers
(`Gene_obj_comparator`, `Gene_validator`, `CDNA_stitcher`,
`Gene_obj_alignment_assembler`, `Annot_update`), adapted to LRAA isoforms:

1. **Loci.** Group reference genes and isoforms by genomic overlap per strand;
   overlaps on the opposite strand are recorded for antisense calls. Loci are
   independent, so a round runs in parallel over them.
2. **Classify each isoform** against the genes it overlaps:
   - **compatible** with a gene's model (same or consistent intron chain):
     incorporated as is, or it extends the UTRs, or it changes the protein;
   - **incompatible:** a new alternative isoform, a gene-structure update, a
     merge of several genes, or a split of one gene;
   - **no gene:** a novel gene, antisense, or intergenic noncoding.
3. **Build and validate the proposed model.** The CDS comes from the reference
   model when the isoform keeps it, otherwise from the isoform's TransDecoder ORF.
   The checks are PASA's: minimum overlap with the original protein, percent
   identity of the protein comparison, length relative to the original
   (complete and partial isoforms have separate thresholds), minimum ORF size,
   minimum fraction of the transcript that is coding, maximum UTR exons, consensus
   splice sites. The thresholds keep PASA's names and defaults
   (`MIN_PERCENT_OVERLAP`, `MIN_PERID_PROT_COMPARE`, ...).
4. **Apply** the accepted changes to give the next `AnnotationSet`.

**PASA's FL vs EST paths.** PASA handles full-length cDNA assemblies (they may
replace or add isoforms) differently from EST assemblies (they are stitched into
existing models to extend them). Proposal: an isoform takes the full-length path
when its TransDecoder ORF is complete (start and stop codons), as PASA's own
TransDecoder plug-in already defined full-length; other isoforms take the
stitching path. LRAA's TSS/polyA support is recorded but not required.

**Single-exon isoforms** follow PASA: they can match or extend existing genes, and
become a novel gene only with a complete ORF above the minimum ORF size.

**Expression.** Every decision and output carries the isoform's expression.
Optional thresholds (minimum TPM, minimum isoform fraction) can keep low-expression
isoforms out of the update; by default LRAA's own filtering decides.

### 5.2 Decisions

PASA's 45 status codes collapse to these actions, each with the PASA status ids it
covers recorded in the docs:

| Action | Changes the annotation | PASA statuses |
|---|---|---|
| `matches_existing` | no | 3, 12 |
| `utr_update` | yes | 4, 13 |
| `protein_update` | yes | 6, 8, 14, 16 |
| `new_isoform` | yes | 9, 24, 25, 27, 33, 44 |
| `merge_genes` | yes | 29, 36 |
| `split_gene` | yes | 40 |
| `replace_unsupported_model` ("stomp") | yes | 26, 31, 32 |
| `novel_gene` | yes | 10 |
| `antisense` | no (reported) | 42, 43 |
| `intergenic_noncoding` | no (reported) | 19 |
| `rejected` | no | 1, 7, 15, 18, 21–23, 28, 30, 34, 35, 41, 45, with the failed check as the reason |
| `deferred` | no, retried next round | 37–39 |

PASA's deprecated statuses (2, 5, 11, 17, 20) are dropped.

### 5.3 Rounds

LRAA's isoforms are fixed. Each round compares every isoform not yet incorporated
with the annotation produced by the previous round. An isoform that could not be
placed earlier may fit once a neighboring change has been made (an extended gene,
a merge, a split); PASA's "delayed incorporation" cases become ordinary retries.

- Round *k* reads `AnnotationSet` *k-1* and writes `rounds/round_k/`: the
  annotation, the decisions table, and a summary.
- Each isoform's state is `incorporated` (round and action), `rejected` (latest
  reason; retried in later rounds) or `pending`.
- Round 1 evaluates every locus. Later rounds revisit only loci where the
  previous round changed or created a gene and isoforms are still unincorporated;
  every other locus would get the same answer again.
- Stops after `--max_rounds` (default 3), or as soon as a round incorporates
  nothing new.
- Rounds remember which isoforms built or confirmed each model. A model built from
  a full-length isoform is not replaced by an incompatible isoform in a later round;
  that isoform is tried as an alternative isoform instead. (Separate Perl PASA runs
  forget this, so a second run can overwrite what the first built.)
- Final outputs: `annotation.updated.gff3` (+ GTF, BED, cDNA/CDS/protein FASTA) and
  `annot_update.decisions.tsv` (every isoform, its final state, the round it
  happened in, the gene(s) affected, the checks, its expression).

**IDs.** Reference genes and transcripts keep their ids; changed models carry a
`pasa_xg_update` attribute and the supporting isoform id. New isoforms are named
`<gene_id>.pxg<N>`, novel genes `PXG_gene_<N>` (transcripts `PXG_gene_<N>.t<M>`). A merge keeps the leftmost
gene's ids (attribute `merged_from`); a split keeps the original ids for its first
product and names the others `PXG_gene_<N>` (attribute `split_from`).

## 6. Gene models from isoforms

Without an annotation (or for novel loci), coding gene models come straight from
the isoforms: genes from LRAA's gene grouping, transcripts from isoforms with a
TransDecoder ORF above the minimum length, CDS mapped to the genome, partial ORFs
labeled. This replaces PASA's `alignment_assembly_to_gene_models.dbi` and its
`pasa_asmbls_to_training_set` TransDecoder route, and gives a training set for
gene finders.

## 7. Alternative splicing

Events between isoforms of the same gene: skipped exons, retained introns,
alternative 5'/3' splice sites, mutually exclusive exons, and alternative first/last
exons, labeled as affecting UTR or protein (PASA's
`classify_alt_splice_as_UTR_or_protein`). Expression from LRAA lets events be
reported with isoform fractions.

The design waits on two things you have: the Python alt-splicing code to examine
or borrow, and the earlier work on alt TSS/polyA analysis and reporting. PASA's
Perl version (`Alternative_splice_comparer.pm`, `classify_alt_splice_isoforms*`,
`splicing_variation_to_splicing_event.dbi`; about 3,600 lines) is the fallback
reference.

## 8. Run directory, caching and resume

Same approach as LRAA's caches and checkpoints:

```
run/
  manifest.json            # PASA-XG version, inputs, settings, tool versions, stage status
  00_prep/  01_align/  02_assemble/  03_isoforms/  04_coding/
  05_annot_update/rounds/round_01/ ...
  06_gene_models/  07_alt_splice/
  outputs/                 # final deliverables
```

- Each stage writes a `stage.done.json`: a digest of its inputs (content hashes of
  small files; size, mtime and a sampled hash for large ones), its settings and
  tool versions. A stage reruns when the digest no longer matches and is skipped
  otherwise, so a run resumes where it stopped and a stale cache is never reused.
- `--from <stage>` forces a stage and everything after it to rerun.
- Expensive reusable artifacts (the minimap2 index, the parsed annotation) can be
  shared through `--cache_dir`, keyed by content digest.
- Tables are Parquet (fast, typed; `pyarrow`). Model sets are pickled for speed,
  with GFF3 always written alongside, so every checkpoint can also be read by
  people and other tools.
- Within a stage, data stays in memory; parallelism is by locus (annotation update,
  alt splicing) or by chunk (TransDecoder).

## 9. Interface and packaging

```bash
pasa-xg run --genome g.fa --reads reads.fq.gz --input_type pacbio_hifi \
            --annotation annot.gff3 --cpu 16 --out run/
pasa-xg run --genome g.fa --bam aligned.bam --input_type ont --out run/
pasa-xg <stage> --out run/          # one stage alone: prep, align, assemble, isoforms,
                                    # coding, annot-update, gene-models, alt-splice, report
pasa-xg status --out run/           # stage states from the manifest
```

- Python ≥ 3.10, installable with pip; the package is `pasa-xg`, the module
  `pasa_xg`.
- External tools: minimap2, samtools, LRAA, TransDecoder; optionally hmmscan,
  BLASTP or DIAMOND. Versions are recorded in the manifest.
- A Docker image with everything installed, and a WDL workflow (Terra/miniwdl),
  as for LRAA.

## 10. What carries over

- **From PyPASA:** the project layout, the CLI and output patterns, and the FASTA
  and GFF writers. Its alignment, validation, clustering and assembly stages are
  replaced by minimap2 and LRAA.
- **Ported from Perl PASA:**

| Perl source | Lines | Ported into |
|---|---|---|
| `cDNA_annotation_comparer.dbi` | 3,511 | annotation update, §5.1 |
| `cDNA_annotation_updater.dbi`, `PASA_UPDATES/*.pm` | 1,100 | applying updates, §5.1 step 4 |
| `Gene_obj.pm` (needed parts), `Gene_obj_comparator.pm`, `Gene_validator.pm` | 6,300 | data model and checks, §4 and §5.1 |
| `CDNA_stitcher.pm`, `Gene_obj_alignment_assembler.pm` | 430 | stitching path, §5.1 |
| `alignment_assembly_to_gene_models.dbi` | 174 | §6 |
| alt-splice scripts and `Alternative_splice_comparer.pm` | 3,600 | §7, alongside your Python code |

- **From LRAA:** branch `pasa_mode_v2` (`--pasa_mode`, `--preassembled_input`,
  chunked-mode fixes). It goes into `devel` after an independent review confirms
  that runs without the new flags are unaffected; PASA-XG then pins the LRAA
  release that contains it.

## 11. Validation

Results are not expected to match Perl PASA. The checks are:

1. **Unit tests per decision type:** small synthetic loci (a reference gene plus
   isoforms built to be a UTR extension, a new isoform, a merge, a split, a novel
   gene, antisense, a failing protein change), as in LRAA's `testing/pasa_mode`
   fixture.
2. **Sample data vs Perl PASA:** run Perl PASA's annotation comparison
   (`-A`, `orig_annotations_sample.gff3`) on the sample data with minimap2 for
   reference. Compare at the decision level and explain the differences; they
   need not be zero.
3. **Recovery benchmark:** degrade a good annotation in known ways (drop
   isoforms, trim UTRs, merge and split genes, delete genes), run PASA-XG, and
   measure how much is restored, with Perl PASA on the same inputs as the
   comparison point. Datasets: human (GENCODE; a few chromosomes) with PacBio
   Kinnex/HiFi and with ONT reads; Arabidopsis (Araport11) with public long
   reads; and Trinity assemblies for `--input_type transcripts`. This doubles as
   a figure for the paper.
4. **End-to-end test:** a tiny dataset through every stage in the Docker image,
   including resume after interruption and `--from`.

## 12. Milestones

| | Milestone | Done when |
|---|---|---|
| M0 | Repo, CLI skeleton, run directory, manifest and stage caching | `pasa-xg run` resumes and skips correctly on stub stages |
| M1 | `prep`, `align`, `assemble`, `isoforms` | sample data and a long-read test set reach the isoform table |
| M2 | `coding` and `gene_models` | TransDecoder ORFs on the genome; gene models GFF3 validated |
| M3 | Annotation update, one round | unit tests for every action; sample-data comparison with Perl PASA reviewed |
| M4 | Multiple rounds, final outputs, `report` | rounds converge on the sample data and the benchmark |
| M5 | Alternative splicing | after reviewing your code and the alt TSS/polyA work |
| M6 | Docker, WDL, benchmark, docs | benchmark results written up |

## 13. Decisions (2026-10-25)

1. Full-length path: a complete TransDecoder ORF; TSS/polyA support not required (§5.1).
2. LRAA `pasa_mode_v2`: independent review first, then merge into `devel` (§10).
3. Single-exon isoforms: novel gene only with a complete ORF above the minimum size (§5.1).
4. IDs: `<gene_id>.pxg<N>`, `PXG_gene_<N>` (§5.3).
5. Rounds: default 3; later rounds revisit only loci changed in the previous round (§5.3).
6. Strand: no option; minimap2 `-ub` for every input, and LRAA infers transcribed
   strands from `ts` on its own (§2).
7. Benchmark: human with PacBio and ONT, Arabidopsis, Trinity assemblies (§11).

Decisions (2026-10-26):

8. LRAA `pasa_mode_v2` merged into `devel` (2a6524d) after review.
9. TransDecoder decides complete vs partial ORFs; no start-codon override. 5'-partial
   isoforms take PASA's non-FL (EST) path.
10. Supported-model protection stays on by default (`--replace_supported_models` for
    PASA's behavior).
11. `MIN_PERCENT_PROT_CODING` stays at 40% (conservative), plus the NMD-target rule:
    no update may yield a stop codon > 50 nt upstream of the last exon-exon junction.
12. PASA's `TRUST_FL_STATUS` on by default: LRAA trims isoform ends, so a compatible
    FL isoform shorter than the model must still count as FL support.

Still open: the alternative-splicing design (§7), pending your code and the alt
TSS/polyA work.
