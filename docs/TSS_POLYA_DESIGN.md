# TSS / PolyA-aware annotation update: design draft

Status: draft for review, 2026-11-02. Nothing here is implemented yet.

## Problem

PASA's update treats transcript ends as unreliable, which made sense for ESTs:

- `merge_compatible` takes the union of the model and the isoform. An end-only difference is
  always a UTR update toward the **outermost** end (statuses 4 / 13). The update only extends;
  it never shortens or adds an isoform.
- An isoform with a shorter but real TSS or polyA site is folded into the longer model, and its
  site is lost.
- An isoform reaching past a real model end replaces it, so the model's own site is lost.
- New isoforms come only from intron-chain or protein differences. Two models differing only at
  their ends are never both kept.
- New isoforms are de-duplicated on exact coordinates only.

PASA-XG already reads LRAA's `TSS` / `PolyA` flags into the isoform table
(`isoforms.py`, also from `--isoforms_gtf`), but `annot_update` never looks at them.

The E5 end calibration shows why the outermost end is the wrong default. On GENCODE-matching
chr22 chains, the longest read end lands within 10 bp of a CAGE peak 24% of the time. The f20
end does so 57% of the time, and each sample's own f20 end 77%.

## What LRAA does (merge_LRAA_GTFs.py, Splice_graph, Simple_path_utils)

- TSS and polyA sites are splice-graph nodes. Sites within `max_dist_between_alt_TSS_sites`
  (50 bp), and the same for polyA, are one site.
- A path that starts on a TSS node or ends on a polyA node is closed at that end. A longer path
  doesn't contain it, so it stays a separate isoform.
- An end with no site is open. A compatible longer path absorbs that isoform.
- Ends near a site are trimmed onto it (`refine_TSS_simple_path`, `refine_PolyA_simple_path`).
- When calling sites from reads, a weak alternative site inside a stronger one is treated as
  degradation and pruned (`max_frac_alt_TSS_from_degradation` 0.20, `min_TSS_iso_fraction`
  0.05). The merge turns these filters off, because there "support" counts input models, not
  reads.

## Proposal

### 1. Sites

- **LRAA:** an isoform end flagged `TSS "True"` or `PolyA "True"`, whether from assembly or
  from `--isoforms_gtf`.
- **Known sites:** `--tss_sites` and `--polya_sites` take BED files with strand, such as FANTOM5
  CAGE peaks, refTSS, PolyASite or PolyA_DB.
- **Window:** `--site_window`, default 50 bp to match LRAA.
- **Supported end:** an isoform end is supported if LRAA called it, or if it lies within the
  window of a known site on the same strand.
- **Anchored model end:** a reference model's end is anchored if it lies within the window of a
  supported isoform end or a known site. Under `--additive`, every reference end is anchored.

### 2. End-only differences

These are isoforms compatible with a model whose intron chain the model already contains or
matches, today's statuses 3, 4, 12 and 13. Each end of the isoform is compared with the model's
end. Within the window, it is the same end.

| isoform end | model end anchored | model end not anchored |
|---|---|---|
| supported, beyond the model end | **new isoform** (alternative TSS / polyA) | UTR update to the isoform's end (as now) |
| supported, inside the model end | **new isoform** (alternative TSS / polyA) | see question 2 below |
| not supported, beyond the model end | model kept (no extension) | UTR update (as now) |
| not supported, inside the model end | absorbed (as now) | absorbed (as now) |

- An alternative-end isoform has the model's CDS with a different UTR. It passes the protein
  check by construction, and its UTR exons still go through `check_utr_exons`.
- Weak alternative ends: an alternative-end isoform must reach `--alt_end_min_fraction`
  (default 0.20) of the TPM of the isoform supporting the model's end. This is LRAA's
  degradation rule, applied to isoforms. With `--isoforms_gtf` there is no TPM, so the test is
  skipped.

### 3. Intron-chain changes, new isoforms and novel genes

- **Taking ends:** when a FL isoform replaces or adds a structure (statuses 6, 8, 9, 10, 44), it
  takes its own ends when they're supported.
- **Unsupported ends:** an unsupported end keeps the model's end when the model's is anchored and
  compatible. Otherwise it keeps the isoform's end, as now.
- **Stitching (non-FL, status 16 / 25):** model ends survive unless the isoform's end at that
  side is supported.

### 4. De-duplication

When a new isoform is added, two transcripts count as the same isoform if they have the same
introns and CDS and both ends fall within the window. The one whose ends are supported is kept,
or the existing one on a tie. This stops end-variants from piling up, which matters for E5,
where merged models multiply with sample count.

### 5. Outputs

Each transcript gets `tss_support` and `polya_support` attributes:

- `lraa`
- `known:<bed name>`
- `reference`
- `none`

Optionally it also gets the matched site coordinate. `annot_update.decisions.tsv` gains new
action types: `alt_tss_isoform`, `alt_polya_isoform` and `end_kept_anchored`.

### 6. Defaults

- **Default mode:** every rule above.
- **`--additive`:** every reference end is anchored, so an existing end never moves. Alternative
  ends become new isoforms, and missing UTRs are still added where the model has none at that
  end (E0's UTR-trim case). A truncated CDS-only end, with no UTR at all, is not treated as
  anchored even under `--additive`.
- **No sites at all:** no LRAA flags and no BEDs means today's behavior exactly, so old results
  stay reproducible.

## Questions for review

1. **Snapping.** Should a supported end move onto the site (the LRAA node coordinate, or a CAGE
   peak's dominant position), or keep the isoform's own coordinate? The draft keeps the isoform
   coordinate and records the site in an attribute.
2. **Shortening.** When the model end is not anchored and a FL isoform has a supported end
   inside it, should the model be trimmed back to the site? PASA never shortens. The draft does
   it in default mode only, never under `--additive`.
3. **Reference ends in default mode.** Should a reference end with no nearby site still count
   as anchored (trusting the curators), or only with evidence? The draft requires evidence in
   default mode.
4. **Weak alternative ends.** Is 0.20 of the TPM supporting the main end a good default for
   keeping an alternative end? LRAA uses 0.20 when calling sites from reads.

## Evaluation

- **E0 chr22 and genome-wide (human PacBio and ONT):** the ends of recovered UTRs and
  alternative isoforms against GENCODE truth, plus the unchanged-model rate.
- **TUSCO:** reference-match status depends on ends within 50 bp, so this is the direct test.
- **E5 chr22 pilot** (11- and 101-sample merges into GENCODE v50): how many alternative-end
  isoforms are kept, and how many of their ends match CAGE / PolyASite against the null
  (200 bp shifted) rate.
- **Arabidopsis E0:** a check outside human. Araport11 ends and PlantAPAdb or TAIR polyA
  sites, if available.
