# Noncoding isoforms and an evidence filter: design draft

Status: agreed 2026-11-04; implemented as `--include_noncoding` (pasa_xg/annot/noncoding.py, see
"Implementation" at the end). The
default behavior stays Perl PASA's. Calibration details: PASA-XG-Paper
`experiments/E5_gencode_xg/results/FILTER_CALIBRATION.md`.

## Problem

PASA's update is coding-only. In the E5 10-sample chr22 pilot (additive mode), 1,257 isoforms were
incorporated into protein-coding genes, but only 16 into noncoding-only GENCODE genes; 824 were
rejected there and 642 were called antisense. The ways a real transcript is turned away:

| status | reason | PASA-XG today |
|---|---|---|
| 1 | protein coding below `MIN_PERCENT_PROT_CODING` (40%) | rejected |
| (demotion) | no ORF: FL isoform treated as non-FL | usually rejected later |
| 30 / 35 | "a matched gene is noncoding" (compare.py:700) | rejected |
| 19 | intergenic, no ORF | reported, never incorporated |
| 18 / 28 etc. | NMD target, stitched model fails validation | rejected |

Two more problems turned up in the pilot:

- **Antisense mislabel.** `classify_antisense` (compare.py:761) re-checks statuses 1, 10, 44, 19, 22,
  27 and 34 on the opposite strand. It relabels an isoform antisense whenever any opposite-strand
  gene overlaps, even when the isoform already matched a same-strand gene. In the pilot, 1,201 of
  the 1,914 "antisense" isoforms share introns with a same-strand gene.
- **Garbage in.** PASA assumes every incoming isoform is a candidate. A per-sample LRAA merge is
  not that: 48% of the pilot's merged isoforms come from one sample. Once noncoding isoforms are
  admitted, the ORF and protein checks no longer screen anything out, so the evidence filter
  below has to come first.

## Evidence (E5, 101 Isopedia samples, chr22 train / chr21 held out)

Label: the intron chain equals a GENCODE v50 chain (a positive). Novel chains are unlabeled, not
negatives, so the "ref recovered" figures measure sensitivity on known isoforms, and "novel kept"
is an upper bound on what a rule lets through.

- All spliced merged isoforms are canonical; LRAA already drops non-canonical introns.
- Recurrence carries most of the signal. The chain-match rate rises from 4.6% for 1 sample to 10%
  for 2, 18% for 5–9 and 37% for 20+ (chr22). Single-sample isoforms stay at 2–5% whatever their
  read count.
- End support helps single-sample isoforms. With neither end within 25 bp of a CAGE peak or a
  PolyASite cluster, 2–7% match; with one or both ends supported, 8–16% match, about the rate of
  2-sample isoforms.
- A-rich genomic sequence past the 3′ end (internal priming) halves PolyASite support (11–18% vs
  33%) but does not change the chain-match rate. It is a reason to distrust a 3′ end, not an
  isoform.
- Mono-exonic isoforms: 867 of 1,251 on chr22 lie inside a reference intron. They match a
  reference single-exon transcript < 1% of the time unless they occur in ≥ 10 samples (20%).

Held-out chr21 (11,846 spliced isoforms; 1,857 reference chains, 8,343 novel):

| rule | ref recovered | novel kept |
|---|---|---|
| n_samples ≥ 2 | 84.9% | 4,801 |
| **n_samples ≥ 2, or ≥ 1 supported end** | **93.1%** | **5,656** |
| n_samples ≥ 3, or (2 and 1 end), or (1 and both ends) | 81.6% | 4,277 |
| n_projects ≥ 2 | 82.6% | 4,599 |
| PU classifier p ≥ 0.2 (all features) | 94.3% | 5,540 |
| PU classifier p ≥ 0.3 (all features) | 88.4% | 4,483 |
| PU classifier p ≥ 0.3 (evidence features only) | 89.7% | 4,951 |
| SQANTI3 ML filter, defaults (chr22 only) | 32% | — |

- The PU classifier is a gradient-boosted model trained on chr22, positives versus unlabeled, with
  the Elkan–Noto correction. On held-out data, ROC AUC is 0.77 with all features and 0.73 with
  evidence only.
- At the same sensitivity, the classifier keeps about 10% fewer novel chains than the simple rules.
  Its top features are n_samples, n_tissues, platform mix and TSS support, the same evidence the
  rules use.
- Dropping the structural features costs little. The model is not mainly learning what annotated
  structures look like.
- SQANTI3's default filter rejects most GENCODE matches. It is trained on reference-match versus
  NNC, and it learns that ends near reference ends are good. It keeps almost all mono-exonic
  isoforms, because those are outside its model.

Conclusion: use a simple, explainable rule as the front end, and keep the classifier score as an
optional annotation and a check on the rule.

## Proposal

### 1. Evidence filter (gates the noncoding pass)

PASA-XG reads per-isoform evidence from GTF transcript attributes when they are present:
`n_samples`, `n_projects`, `TSS`, `PolyA`, `percA_down`. A helper writes them onto the merged GTF
from the merge tracking file (the paper's `evidence.py`, moved into PASA-XG as
`pasa_xg evidence`).

- The thresholds are relative to cohort size: `n_samples` counts are given as `max(k, f × cohort)`.
  At 101 samples, f is small and k sets the floor.
- The filter applies to the noncoding pass only (decided 2026-11-04). The coding path is PASA's,
  unfiltered. Among isoforms that pass PASA's coding checks, the rule did not separate known from
  new: in the 10-sample pilot it would drop 21% of GENCODE matches and 18% of new isoforms.
- An isoform that the noncoding pass turns away for evidence keeps its coding-path status, with
  `filtered: <reason>` added to its comment.

| noncoding-pass class | default admission |
|---|---|
| ORF consistent with the gene's product, not NMD (§5) | n_samples ≥ 2, or ≥ 1 supported end |
| other spliced | n_samples ≥ 2 and ≥ 1 supported end, or n_samples ≥ 3 |
| mono-exonic | n_samples ≥ 3, not inside a reference intron, and both ends supported |

The mono-exonic rule is provisional and needs research. In a merged catalog, a mono-exonic
isoform's n_samples counts samples whose isoform lies inside it, usually one near-full-length
isoform plus shorter ones nested in it. A count of samples covering ≥ 80% of the span is the
candidate replacement. Full coding-path filtering stays a possible opt-in (`--filter_coding`),
to test on the 101-sample catalog.

Without these attributes, as with a single LRAA run, the filter is off, and a warning says so.

### 2. `--include_noncoding`

After the coding rounds, a noncoding pass takes isoforms that were turned away only for coding
reasons (statuses 1, 19, 30/35 "matched gene is noncoding", the no-ORF demotion, and NMD failures)
and that pass the noncoding admission rule.

- **Known gene:** compare by structure only, against all of the gene's transcripts, noncoding
  ones included. Reuse `merge_compatible`: if the isoform is compatible with an existing noncoding
  transcript, it is a UTR/end update of that transcript; otherwise it is a new isoform. Coding
  models are never changed by a noncoding isoform.
- **Novel locus:** spliced isoforms only, as a new gene.
- **Labels on transcripts:**
  - `transcript_type "pasa_noncoding"`;
  - `pasa_xg_coding "none" | "low" | "nmd"`, with `pasa_xg_orf_len` when an ORF exists.
- **Label on genes:** `pasa_xg_gene_class "likely_noncoding"` when a gene has no coding transcript,
  or `"coding"` otherwise. GENCODE's own `gene_type` is kept.
- A coding/noncoding classifier (e.g. CPAT or RNAsamba) can later add a `coding_potential`
  attribute; it does not gate admission.

### 3. Antisense relabel fix

`classify_antisense` relabels only isoforms that matched no same-strand gene: intergenic (19),
novel gene (10 / 22) and novel alternative isoforms (44). Statuses 1, 27 and 34, where a
same-strand gene was matched, keep their label. This changes reported counts only, not the
annotation, so it can go in default mode too. Perl parity tests compare annotations, not this
label.

### 4. Parity

Without `--include_noncoding` and `--isoform_filter`, the outputs are identical to today's. The
parity tests run as now.

### 5. Coding-like isoforms in coding genes: a lower bar

The user's suggestion (2026-11-04): a noncoding-pass isoform in a known coding gene can be held to
lower thresholds when it is not an NMD target and its ORF is consistent with the gene's expected
product. "Consistent" here means the ORF's protein aligns to one of the gene's reference proteins
under the existing `compare_proteins` test, or is a contiguous part of one: same frame, shared
start or stop.

Such an isoform is labeled `pasa_xg_coding "low"` with its ORF kept as the CDS, and is admitted
under the coding-path rule (n_samples ≥ 2, or ≥ 1 supported end) rather than the stricter
noncoding rule. It typically has a coding fraction under 40%, i.e. a long UTR or retained UTR
intron. NMD targets and isoforms whose ORF is unrelated to the gene's protein stay at the
noncoding bar.

## Decisions

Agreed 2026-11-04:

- Admission thresholds as in §1 (the recommended rule). Keep the classifier results on record and
  revisit later.
- The noncoding second pass (§2), including noncoding isoforms of known coding genes.
- §5: a lower bar for non-NMD isoforms whose ORF fits the gene's product.

Also agreed 2026-11-04:

- The evidence filter gates the noncoding pass only; the coding path is unchanged.
- The antisense relabel fix (§3) applies in default mode.
- The never-shorten guard applies in default mode when TSS / polyA sites are given: a compatible FL
  isoform whose merged model fails the protein comparison no longer cuts the model back to the
  isoform. Without sites, PASA's behavior is kept.

Open: the mono-exonic rule (research pending).

## Implementation (2026-11-04)

`--include_noncoding` runs `annot/noncoding.py` after the coding rounds, on the final annotation,
over every eligible isoform that the coding rounds did not incorporate (deferred ones aside). Its
decisions are one more round in `annot_update.decisions.tsv`. Statuses: 101 added to a known gene,
102 novel noncoding gene, 103 not enough evidence, 104 assembles with an existing transcript,
105 two or more genes, 106 mono-exonic with no gene. The comment keeps the coding-round outcome.

- **Gene assignment** (`assign_genes`). Among same-strand genes sharing an intron with the
  isoform, take the gene whose best transcript's intron set is most similar (Jaccard); ties are
  status 105. This resolves GENCODE readthrough genes (PRR5 / PRR5-ARHGAP8) and overlapping
  loci. With no shared intron, use exon overlap, ignoring genes whose transcripts are all
  < 300 bp (embedded miRNA / snoRNA genes). On the 10-sample pilot, multi-gene skips went from
  1,603 (exon overlap only) to 135.
- **Redundancy.** An isoform that `merge_compatible`s with any transcript of its gene, coding or
  not, adds nothing (104). No end updates of noncoding transcripts in this version.
- **Labels.**
  - `pasa_xg_coding`: `low` when the isoform's ORF (its CDS from the coding stage) is a substring
    of a gene protein or passes `compare_proteins` (the CDS is kept); `nmd` when the stop codon is
    more than `NMD_STOP_TO_JUNCTION` upstream of the last junction; else `none`.
  - `transcript_type`: `pasa_low_coding` or `pasa_noncoding`.
  - Every gene gets `pasa_xg_gene_class` `coding` / `likely_noncoding` (on GTF transcript lines
    too). Novel genes are `PXG_gene_N`, with `pasa_xg_novel "noncoding"` and `antisense_to` when
    they overlap a gene on the other strand.
- **Evidence.**
  - `n_samples` comes from the isoform GTF attribute. When it is missing, it is taken as 1, so
    only end support can admit, rather than the filter being off.
  - Supported ends: the isoform's TSS / PolyA flag, or a known site within `--site_window`
    (default 50 bp; the calibration used 25 bp).
  - For merged catalogs, `tag_ends.py --tracking` (PASA-XG-Paper) writes `n_samples`.
- **Pilot** (E5 random-10, chr22, additive): 1,108 isoforms added to known genes (625 low, 244
  NMD, 239 none), 68 in novel noncoding genes. Antisense calls in round 1 fell from 1,861 to 663
  with the relabel fix.
