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

Final annotations (PASA-XG after 3 rounds vs Perl after its 2 rounds, as in PASA's
sample pipeline): genes 439 vs 456; original genes changed 139 vs 131 (118 in common);
new or changed structures 375 vs 378 (246 identical). Perl's second round mostly
re-confirms its first (3 new structures); PASA-XG's round 2 incorporates 21 more.

Disagreements inspected so far follow from different assemblies (LRAA vs PASA's
assembler) feeding the same validation rules, not from a rule difference.

Perl reference run: `PASA_port/perl_ref_annot/run_annot.sh` (needed DB_File and FASTA36
added to `tools/perl_env`).

## PacBio test (LRAA testing/single_contig: 2,507 reads, GENCODE subset of 50 genes)

`--bam` input, `--input_type pacbio_hifi`: 115 isoforms; round 1 incorporates 69
(36 protein updates, 14 new isoforms, 13 matches, 6 UTR updates), round 2 none.
This run exposed 5'-incomplete reference CDSs (GENCODE `cds_start_NF`) being translated
from the wrong frame; the first CDS phase is now honored everywhere (0 of 336 reference
CDSs flagged, was 89).

## Not done yet

- Alternative splicing (waiting on the user's Python code and alt TSS/polyA work).
- Docker image and WDL: drafts in `docker/` and `WDL/`, not built or run.
- Benchmarks (human with PacBio/ONT, Arabidopsis, Trinity).
- Repo not yet on GitHub (`PASApipeline/PASA-XG` to be created).

## Recovery benchmark: human chr22, PacBio Kinnex PBMC reads (2026-10-26)

Truth: GENCODE47 chr22 protein-coding genes (428). Degraded input (scripts/bench_degrade.py,
seed 1): canonical transcript only; 48 genes deleted, 87 cut back to their CDS, 21 split
in two. Scored by scripts/bench_score.py ("supported" = an LRAA isoform has the truth
structure, i.e. what the evidence allows). Data under bench/ (not in git).

Whole chr22 (1.29M reads, `--bam`, 8 CPUs): 29 min, 1.35 GB peak; LRAA is 27 of the 29 min.

| | PASA behavior (`--replace_supported_models`) | default | default + `--min_fl_orf_aa 100` |
|---|---|---|---|
| deleted genes recovered (25 supported) | 21 | 21 | 25 |
| UTR-trimmed models regaining UTR (41 supported) | 45 | 44 | 48 |
| split genes rejoined with truth structure (8 supported) | 7 | 8 | 9 |
| alternative isoforms restored (277 supported) | 96 | 102 | 110 |
| correct canonical models changed (of 272) | 97 | 13 | 29 |
| added multi-exon structures found in GENCODE | 42% | 34% | 31% |

Region chr22:42-45 Mb (83,553 reads, same reads as FASTA for both pipelines):

| | Perl PASA | Perl, MIN_FL_ORF_SIZE=100 | PASA-XG | PASA-XG, --min_fl_orf_aa 100 |
|---|---|---|---|---|
| wall time / peak memory | 77 min / 4.3 GB | | 6 min / 1.4 GB | |
| deleted genes recovered (3 supported) | 0 | 0 | 2 | 3 |
| UTR-trimmed models regaining UTR (4 supported) | 0 | 0 | 4 | 4 |
| alternative isoforms restored (16 supported) | 1 | 4 | 5 | 5 |
| added multi-exon structures (in GENCODE) | 1 (0) | 27 (4) | 60 (13) | 86 (17) |

Perl PASA (all reads marked full-length) chains long reads into very large assemblies
(up to 839 kb) and rejects 820 of 899 for low coding fraction; PASA's assembler is not
built for long reads, which is the case for PASA-XG.

Fixes found by the benchmark: TransDecoder 5'-partial ORFs on full-length isoforms (now
started at the first in-frame ATG); supported reference models no longer replaced
(`--replace_supported_models` restores PASA's behavior).

### Revision (2026-10-26): TransDecoder's ORF calls, NMD rule, TRUST_FL_STATUS

- TransDecoder decides complete vs partial; the first-ATG override is gone, so a
  5'-partial ORF takes PASA's non-FL (EST) path.
- `REJECT_NMD_TARGETS` (default on): no update may produce a model whose stop codon lies
  more than 50 nt upstream of its last exon-exon junction.
- `TRUST_FL_STATUS` (PASA's `--TRUST_FL_STATUS`, default on): a compatible FL isoform
  that the model extends past stays FL. PASA otherwise treats it as an EST, so it no
  longer shields the model from replacement (status 8); with LRAA's trimmed ends that
  was the cause of most damaged models (58 of 91 cases; 21 more: supporter < 40% coding).
- Percent coding (CDS / transcript), chr22: GENCODE canonical 32% below 40% (11% below
  20%); LRAA complete-ORF isoforms 41% below 40% (10% below 20%); NMD-biotype 69% below 40%.

chr22 recovery (same LRAA assembly for all; "of" = supported by the reads):

| | PASA rules | p40, PASA rules, no refine | p40 +NMD +trust_FL | same, no protection | p30 | p20 | p0 | of |
|---|---|---|---|---|---|---|---|---|
| deleted genes recovered | 21 | 16 | 17 | 16 | 19 | 21 | 21 | 25 |
| UTRs restored | 45 | 46 | 45 | 46 | 45 | 47 | 47 | 41 |
| split genes rejoined (true chain) | 7 | 9 | 10 | 9 | 11 | 11 | 11 | 8 |
| alt isoforms restored | 96 | 106 | 107 | 102 | 112 | 113 | 114 | 277 |
| correct models, chain changed (of 272) | 97 | 82 | 8 | 38 | 11 | 11 | 11 | |
| added/changed multi-exon (in GENCODE) | 742 (309) | 754 (284) | 687 (267) | 673 (270) | 762 (296) | 816 (314) | 851 (324) | |

p30/p20/p0 columns: NMD + trust_FL + protection on. ("PASA rules" = the earlier run with
first-ATG starts.)
