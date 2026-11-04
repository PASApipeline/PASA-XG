# Noncoding isoforms and an evidence filter: design draft

Status: agreed in part, 2026-11-04 (see "Decisions" at the end). Not yet implemented. The
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

### 1. Evidence filter (front end, `--isoform_filter`)

PASA-XG reads per-isoform evidence from GTF transcript attributes when they are present:
`n_samples`, `n_projects`, `TSS`, `PolyA`, `percA_down`. A helper writes them onto the merged GTF
from the merge tracking file (the paper's `evidence.py`, moved into PASA-XG as
`pasa_xg evidence`).

- The thresholds are relative to cohort size: `n_samples` counts are given as `max(k, f × cohort)`.
  At 101 samples, f is small and k sets the floor.
- An isoform that fails the filter is reported with action `filtered` and is not evaluated.

| class | default admission |
|---|---|
| spliced, coding path | n_samples ≥ 2, or ≥ 1 supported end |
| spliced, noncoding path (§2) | n_samples ≥ 2 and ≥ 1 supported end, or n_samples ≥ 3 |
| mono-exonic | n_samples ≥ 3, not inside a reference intron, and both ends supported |

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

Open:

1. Does the evidence filter gate the coding path by default once evidence attributes are present,
   or is it opt-in? Recommendation: on by default when the attributes are present (the coding path
   lets single-sample junk through too), and off when they are absent, with a warning.
2. Antisense relabel fix in default mode, given that it changes no annotation? Recommendation: yes.
3. The never-shorten guard (§2 of TSS_POLYA_DESIGN.md) in default mode with sites: should a
   compatible FL isoform whose merged model fails the protein comparison stop cutting the model
   back to the isoform (now done for `--additive` only, 5d17023)? Recommendation: yes when sites
   are given, since the end rules promise models are not shortened; Perl parity without sites.
