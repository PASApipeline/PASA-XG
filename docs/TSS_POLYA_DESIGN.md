# TSS / PolyA-aware annotation update: design draft

Status: design agreed 2026-11-02 (the user's rule on alternative ends, and the four decisions
at the end). Implementation in progress.

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

The user's rule (2026-11-02): **don't annotate a separate isoform for every alternative TSS or
polyA.** Add one only when the current end isn't supported, or when the alternative end changes
the protein, so that the ORF is shorter or longer than the annotated one.

Models are never shortened (decision 2).

| isoform end | protein | model end anchored | model end not anchored |
|---|---|---|---|
| supported, beyond the model end | same CDS | **no change**; the site is recorded on the model | **end update**: extend to the isoform's end |
| supported, inside the model end | same CDS | **no change**; the site is recorded | **no change**; the site is recorded |
| supported, beyond the model end | different CDS | **new isoform** (alternative-end protein) | **end update**, as a protein update (validated as now) |
| supported, inside the model end | different CDS | **new isoform** (alternative-end protein) | **new isoform** (alternative-end protein) |
| not supported, beyond the model end | same CDS | model kept, no extension | UTR update (as now) |
| not supported, inside the model end | any | absorbed (as now) | absorbed (as now) |

When several isoforms support different sites beyond an unanchored model end, the most-expressed
one sets the new end (decision 4).

- **Different CDS:** this happens, for example, when an alternative TSS falls downstream of the
  annotated start codon, so translation starts at a later ATG. It also happens when an
  alternative polyA site cuts into the CDS, or adds an exon that carries a different stop. The
  new isoform goes through the usual protein checks: percent coding, NMD and
  `compare_proteins`.
- **Recording sites:** an alternative site that doesn't produce an isoform isn't thrown away.
  The model records it as `alt_tss_sites` / `alt_polya_sites`, each with its support source, so
  the evidence stays in the annotation for later analysis. Alternative TSS and polyA analysis
  is a downstream utility, deferred.
- **Expression threshold:** none is needed in this design, because alternative ends no longer
  add isoforms on their own.

### 3. Intron-chain changes, new isoforms and novel genes

- **Taking ends:** when a FL isoform replaces or adds a structure (statuses 6, 8, 9, 10, 44), it
  takes its own ends when they're supported.
- **Unsupported ends:** an unsupported end keeps the model's end when the model's is anchored and
  compatible. Otherwise it keeps the isoform's end, as now.
- **Stitching (non-FL, status 16 / 25):** model ends survive unless the isoform's end at that
  side is supported.

### 4. De-duplication

When a new isoform is added, two transcripts with the same introns and CDS count as the same
isoform, whatever their ends. The one whose ends are supported is kept, or the existing one on
a tie, and the other's sites go onto it as `alt_*_sites`. Section 2 keeps a separate isoform
only when the protein differs. This also stops end-variants from piling up in E5, where merged
models multiply with sample count.

### 5. Outputs

Each transcript gets `tss_support` and `polya_support` attributes, plus the `alt_tss_sites` and
`alt_polya_sites` its isoforms supported. The support attributes take these values:

- `lraa`
- `known:<bed name>`
- `reference`
- `none`

Optionally it also gets the matched site coordinate. `annot_update.decisions.tsv` gains new
action types: `alt_end_protein_isoform`, `end_update` and `alt_site_recorded`.

### 6. Defaults

- **Default mode:** every rule above.
- **`--additive`:** every reference end is anchored, so an existing end never moves. An
  alternative end adds an isoform only when its protein differs. Missing UTRs are still added
  where the model has none at that end (E0's UTR-trim case). A truncated CDS-only end, with no UTR at all, is not treated as
  anchored even under `--additive`.
- **No sites at all:** no LRAA flags and no BEDs means today's behavior exactly, so old results
  stay reproducible.

## Decisions (user, 2026-11-02)

1. **Snapping:** keep the isoform's own end coordinate, and record the matched site in an
   attribute.
2. **Shortening:** never. A model end is only ever extended. A supported site inside the model
   is recorded, or becomes a new isoform when the protein differs.
3. **Reference ends in default mode:** anchored only with evidence, meaning a known site or a
   supported isoform end within the window. `--additive` anchors every reference end.
4. **Which isoform sets the end:** the most-expressed supported isoform. With `--isoforms_gtf`
   there is no TPM, so it is the one with the most source support (merge tracking count), else
   the outermost supported end.

## Evaluation

- **E0 chr22 and genome-wide (human PacBio and ONT):** the ends of recovered UTRs and
  alternative isoforms against GENCODE truth, plus the unchanged-model rate.
- **TUSCO:** reference-match status depends on ends within 50 bp, so this is the direct test.
- **E5 chr22 pilot** (11- and 101-sample merges into GENCODE v50): how many alternative-end
  isoforms are kept, and how many of their ends match CAGE / PolyASite against the null
  (200 bp shifted) rate.
- **Arabidopsis E0:** a check outside human. Araport11 ends and PlantAPAdb or TAIR polyA
  sites, if available.
