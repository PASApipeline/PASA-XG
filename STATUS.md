# PASA-XG status

Last updated 2026-10-25 (overnight build).

## Working

| Stage | What | Checked by |
|---|---|---|
| prep | genome faidx + minimap2 index; annotation (GTF/GFF3) parsed and CDS-checked | unit tests; sample data (281 genes, 2 with CDS problems) |
| align | minimap2 per `--input_type`, or a given BAM | sample data |
| assemble | LRAA (`pasa_mode_v2` @ 2a6524d via `scripts/lraa_dev.sh`), ref-guided with an annotation | sample data |
| isoforms | LRAA GTF + expression into one table, isoform FASTA | sample data |
| coding | TransDecoder LongOrfs/Predict (`-S`, single best), ORFs mapped to the genome | all 785 genome-mapped proteins equal TransDecoder's |
| annot_update | PASA's comparer/updater over LRAA isoforms, multi-round | unit tests per decision type; sample data vs Perl PASA (below) |
| gene_models | coding models straight from the isoforms | sample data |
| report | deliverables in `outputs/`, `summary.tsv` | sample data |

Resume: each stage's checkpoint is keyed on its inputs' content, its settings and its
own source code; a rerun stage starts from an empty directory (TransDecoder and LRAA keep
their own checkpoints, which would otherwise serve stale results).

## Sample data (PASApipeline sample_data, transcripts, orig_annotations_sample.gff3)

Whole pipeline: about 2.5 minutes on 4 CPUs (LRAA ~70 s, TransDecoder ~10-90 s).

Annotation update rounds: round 1 evaluates 1006 isoforms and incorporates 602;
round 2 revisits 148 and incorporates 21; round 3 revisits 7 and adds nothing.

Round 1 against Perl PASA's round 1 (same reads, minimap2; Perl uses its own assembler):

| | Perl | PASA-XG |
|---|---|---|
| genes after update | 456 | 440 |
| original genes changed | 131 | 138 (118 in common) |
| isoforms added to original genes | 29 genes | 51 genes |
| new or changed structures (introns + CDS) | 375 | 375 (245 identical) |
| novel genes | 179 | 161 (121 overlap a Perl novel gene) |

Disagreements inspected so far follow from different assemblies (LRAA vs PASA's
assembler) feeding the same validation rules, not from a rule difference.

## Not done yet

- Alternative splicing (waiting on the user's Python code and alt TSS/polyA work).
- Docker image and WDL: drafts in `docker/` and `WDL/`, not built or run.
- Benchmarks (human with PacBio/ONT, Arabidopsis, Trinity).
- Repo not yet on GitHub (`PASApipeline/PASA-XG` to be created).
