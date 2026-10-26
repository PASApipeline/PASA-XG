# Demonstrating PASA-XG as a genome annotation tool: evaluation plan

Draft, 2026-10-26; decisions of the same day are in section 9. Goal: show that PASA-XG keeps a genome annotation current. It turns a
draft annotation into an accurate one, and keeps improving a mature one, as transcript
evidence accumulates: ESTs and full-length cDNAs, short-read assemblies, then long reads.
Alternative-splicing analysis is out of scope for now.

## 1. Claims to support

| | Claim | Main evidence |
|---|---|---|
| C1 | Raises a draft annotation (ab initio, projected, or BRAKER-style) toward reference quality | E1 |
| C2 | Conservative: correct models are left alone; changes are mostly right | E0, E4 |
| C3 | Cross-generational: one engine for ESTs/cDNAs, Trinity/StringTie assemblies, long reads, and their combination | E3 |
| C4 | Maintenance at scale: adding tissues/samples keeps improving the annotation; sequential updates converge to an all-at-once update | E2 |
| C5 | On a reference annotation (GENCODE), population-scale long reads yield a measurably expanded, validated proteome | E5 |
| C6 | Practical: fast, resumable, provenance for every change; complements rather than replaces gene predictors | all; run-time tables |

The positioning that falls out of the literature (section 2): gene predictors build
annotations from scratch, and transcript catalogs list isoforms without a coding-region
model. PASA-XG is the **update layer** between them: it brings any existing annotation
up to date with transcript evidence, adds CDS-consistent models, and refuses changes
that would damage supported models or yield NMD targets.

## 2. Landscape (what exists, what was measured)

**Annotation pipelines**

| Tool | Approach | Transcript evidence | Updates an existing annotation? |
|---|---|---|---|
| MAKER2 ([Holt & Yandell 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3280279/)) | ab initio + evidence alignment, AED-scored | ESTs, RNA-seq assemblies | yes ("can update legacy annotations") |
| BRAKER3 ([Gabriel et al. 2024, Genome Res](https://genome.cshlp.org/content/34/5/769)) | GeneMark-ETP + AUGUSTUS + TSEBRA | short-read RNA-seq + proteins; Iso-Seq "not extensively tested" ([route map](https://arxiv.org/pdf/2403.19416)) | no (de novo) |
| Funannotate ([docs](https://funannotate.readthedocs.io/en/latest/update.html)) | EVM consensus, then **PASA** `update` for UTRs/isoforms | RNA-seq via Trinity + PASA | yes, through PASA |
| EviAnn ([Zimin et al. 2025, bioRxiv](https://www.biorxiv.org/content/10.1101/2025.05.07.652745.full.pdf)) | evidence-only: StringTie2 transcripts + miniprot proteins | short reads / transcripts | no |
| EASEL ([UConn](https://genomics.institute.uconn.edu/biodiversity/software)) | StringTie2/PsiCLASS + TransDecoder + AUGUSTUS, random-forest filter | short reads | no |
| Tiberius ([Gabriel et al. 2024](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11645249/)) | deep learning + HMM, ab initio (mammals) | none | no |
| Helixer ([Holst et al. 2025, Nat Methods](https://www.nature.com/nature-index/article/10.1038/s41592-025-02939-1)) | deep learning + HMM, ab initio (plants, fungi, animals) | none | no |
| Liftoff / LiftOn / TOGA | projection from a related annotated genome | none | projection only |
| Ensembl genebuild | protein, RNA-seq, long-read and projection pipelines | includes long reads | rebuilds |

**Benchmarks and their metrics**

- *Building better genome annotations across the tree of life* ([Genome Res 2025](https://genome.cshlp.org/content/35/5/1261)):
  12 methods, 21 species. TOGA, BRAKER3 and StringTie were the top performers. Metrics:
  BUSCO recovery, CDS length and false-positive rate. RNA-seq "substantially improves"
  annotation when whole-genome alignment isn't feasible.
- BRAKER3 paper: gene/transcript sensitivity and precision against reference annotations in 11
  species; MAKER2 and Funannotate as comparators.
- Tiberius: gene-level F1 on human of 62%, against 21% for the next-best ab initio tool;
  matches BRAKER3, which uses RNA-seq and proteins.
- *Evaluation of strategies for evidence-driven genome annotation using long-read RNA-seq*
  ([Paniagua et al. 2025, Genome Res](https://pmc.ncbi.nlm.nih.gov/articles/PMC12047274/)):
  - Data: LRGASP WTC11 chr19 (against GENCODE) and Florida manatee.
  - Curated PacBio transcript models combined with AUGUSTUS beat ab initio and short-read
    evidence; raw ONT reads did worst.
  - Gains plateau at about half the PacBio reads.
  - Manatee BUSCO: 30.7% ab initio vs 87.9% combined.
  - Against NCBI's manatee annotation, loci agree about 90% but transcripts only about 35%.
  - No method updates an existing annotation; PASA isn't used.
- *Welcome to the big leaves* ([Vuruputoor et al. 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10439824/)):
  - BRAKER and MAKER on five plant genomes.
  - Transcripts alone are insufficient; combining evidence with ab initio is recommended,
    and short plus long reads together help.
  - Protein evidence adds false positives.
  - Metrics: structure statistics, reciprocal similarity search, and mono- vs multi-exonic
    gene counts.
- LRGASP ([Nat Methods 2024](https://link.springer.com/10.1038/s41592-024-02298-3)): for
  transcript identification, read length and accuracy matter more than depth.
  Reference-guided tools do best on well-annotated genomes. Tools agree only moderately,
  and least on poorly studied genes.
- Quality without a truth set: BUSCO/compleasm, and OMArk
  ([Nat Biotech 2024](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11738984/)), which
  measures completeness plus lineage consistency, flagging spurious and fragmented proteins.

**Reference-scale catalogs and protein evidence**

- GENCODE 2025 ([NAR](https://pmc.ncbi.nlm.nih.gov/articles/PMC11701607/)): long reads are
  used to fix and add models. The CLS project
  ([2024](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11565817/)) roughly doubled lncRNA
  genes in v47. Protein-coding curation is still largely manual, aided by proteomics and
  Ribo-seq.
- CHESS 3 ([Varabyou et al. 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10614308/)):
  a GTEx-based catalog of 19,839 protein-coding genes, with 14,863 protein-coding
  transcripts in no other catalog. Filtered by expression, conservation and structure.
- Structure-guided isoform selection ([Sommer et al. 2022, eLife](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9812405/)):
  ColabFold pLDDT on more than 230,000 isoforms; a well-folded isoform argues it is
  functional. This is the precedent for the AlphaFold idea.
- Long-read proteogenomics ([Miller et al. 2022, Genome Biol](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8892804/)):
  long-read protein databases nearly double the protein isoforms detected by mass spec.
- Long-read collections:
  - ENCODE4 ([Reese et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10245583)): PacBio,
    264 libraries, 81 human/mouse samples; at least one full-length transcript for 87.7%
    of human protein-coding genes.
  - GTEx ONT ([Glinos et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC10888753)):
    about 15 tissues; about 70,000 novel transcripts, of which 2,575 have mass-spec support.
- Isopedia ([Zheng et al. 2026, bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.23.713667v1)):
  an index of 1,007 human long-read samples (603 PacBio, 404 ONT, 37 tissues), keyed on
  splice-junction chains. **It is queried with a GTF; it does not export its isoform set**
  ([README](https://github.com/zhengxinchang/isopedia)). The index is 305 GB (107 GB
  compressed) and needs 16 GB of RAM to query.
- Arabidopsis: AtRTD3 ([Zhang et al. 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9264592/)),
  169k transcripts (79% from Iso-seq) and 1,541 genes not in Araport11. RAFL full-length
  cDNAs ([Seki et al. 2002](https://bioweb.supagro.inrae.fr/ESTHER/paper/Seki_2002_Science_296_141))
  and Sanger ESTs. Araport11 itself used a PASA update step.
- LRAA ([preprint](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12919103/)).

**Gaps PASA-XG can fill**

1. No published pipeline takes an existing annotation and updates it with long reads under
   explicit safeguards. Funannotate and MAKER update through PASA or MAKER's own logic,
   built for ESTs and short reads.
2. BRAKER3's long-read support is experimental. Paniagua et al. build anew rather than
   update, and their transcript-level agreement with NCBI was low (about 35%).
3. Catalogs (CHESS, ENCODE4, GTEx-ONT, Isopedia) report transcripts. Turning them into
   CDS-consistent gene models (ORF, start/stop, NMD, protein comparison against the
   current model) is exactly the PASA step.
4. Benchmarks rarely measure damage to correct models or behavior over successive
   updates. Both are central to annotation maintenance.

## 3. Evaluation principles

- **No circularity.** The truth set is never an input or LRAA guide. When the reference
  guides LRAA (E4, E5), success is measured by later curation or orthogonal evidence, not
  by agreement with the reference.
- **Several truth views**, because references are incomplete:
  - (a) exact match to the reference (gffcompare: base/exon/intron/intron-chain/transcript/gene);
  - (b) high-confidence subsets (MANE Select for human; models with full-length cDNA
    support for Arabidopsis);
  - (c) precision against the union of independent catalogs (GENCODE ∪ RefSeq ∪ CHESS 3;
    Araport11 ∪ AtRTD3);
  - (d) reference-free: BUSCO/compleasm, OMArk, CDS length distribution, mono- vs
    multi-exon ratio.
- **Protein level as well as structure:**
  - fraction of reference proteins reproduced exactly;
  - proteins changed vs the start;
  - N-terminal extensions and truncations;
  - Pfam domain gain or loss.
- **Damage:** correct models changed or removed. Keep the chr22 measure, applied genome-wide.
- **Matched evidence:** every comparator gets the same reads and transcripts, with its
  recommended settings.
- **Provenance:** every changed model records which isoform(s), which input set and which
  round changed it. This makes the tissue and generation analyses straightforward.

## 4. Experiments

### E0. Controlled damage and recovery, genome-wide

The chr22 design applied to whole genomes. Truth is the reference itself. Genes are
deleted, UTRs trimmed, genes split and isoforms dropped; then we measure recovery and
damage to the untouched models.

- Species: human (GENCODE MANE subset) and Arabidopsis (Araport11).
- Evidence: one PacBio dataset and one ONT dataset per species.
- Comparators: Perl PASA on the same LRAA isoforms (as done on chr22); Funannotate update.
- Output: recovery vs supported table; damage rate; run time.
- Compute: Arabidopsis whole-genome locally; human per chromosome locally, or all
  chromosomes on Terra.

### E1. From draft to current: the update layer (C1, C6)

Start from annotations a new genome would have, then update each with the same long-read
evidence.

- Starting annotations:
  - ab initio: AUGUSTUS (species model), Helixer, and Tiberius for human;
  - projection: LiftOn/TOGA from a relative (human from mouse; Arabidopsis from
    *A. lyrata* or *Capsella*);
  - full pipeline: BRAKER3 on short-read RNA-seq + OrthoDB proteins.
- Evidence: long reads from one tissue (then E2 for more).
- Comparators, all given the same evidence:
  - BRAKER3 with Iso-Seq (route-map protocol);
  - EviAnn;
  - evidence only: LRAA isoforms + TransDecoder;
  - StringTie2 + TransDecoder;
  - Perl PASA on LRAA isoforms.
- Key figure: for each starting annotation, accuracy before and after PASA-XG, next to
  the comparators. The expected message is that PASA-XG lifts every starting point and
  combines with BRAKER3 rather than competing with it.
- Truth: GENCODE/MANE; Araport11 (+ AtRTD3 for precision). Metrics as in section 3.

### E2. Tissues, merging and incremental maintenance (C4)

- Data:
  - human: ENCODE4 PacBio tissues and GTEx ONT tissues, so both platforms;
  - Arabidopsis: the AtRTD3 Iso-seq tissue and condition libraries.
- Design:
  1. Each tissue alone: LRAA, then PASA-XG.
  2. Cumulative: add tissues one at a time, either (a) LRAA per tissue then
     `merge_LRAA_GTFs.py`, then one update, or (b) sequential PASA-XG updates, each using
     the previous output.
  3. All together.
- Measures:
  - saturation curves (genes and isoforms corrected or added vs number of tissues);
  - tissue-restricted genes recovered only by their tissue;
  - (a) vs (b) concordance, i.e. order independence. That shows the annotation can be
    maintained as new data arrives.

### E3. Cross-generation (C3)

- Arabidopsis is the natural showcase: PASA was first built on its ESTs and full-length
  cDNAs, and every evidence generation exists for it.
  - Gen 1: Sanger ESTs + RAFL full-length cDNAs (`--input_type transcripts`; chaining on
    for ESTs).
  - Gen 2: Illumina RNA-seq → Trinity (genome-guided) and StringTie.
  - Gen 3: Iso-seq / ONT → LRAA.
- Same starting annotation as E1 (ab initio). Each generation alone, then layered in
  chronological order (Gen 1 → +2 → +3), each layer a PASA-XG update of the previous output.
- Measures: what each generation adds, mainly complete CDS, UTRs and isoforms. Overlap
  among the generations' contributions. Long reads should dominate isoforms; cDNAs still
  add genes that are poorly expressed in the sampled tissues.
- Human version, on a few chromosomes: dbEST + MGC full-length cDNAs, then GTEx/ENCODE
  Illumina, then long reads.

### E4. Retrospective "time travel": does PASA-XG anticipate curation? (C2)

- Update an **older** reference with modern long reads, then score each change by whether
  curators made the same change later:
  - human: a GENCODE release from about 2018, scored against v49;
  - Arabidopsis: TAIR10, scored against Araport11 and its later updates.
- Measures:
  - the fraction of PASA-XG changes later adopted by curators, which estimates precision
    independently of the evaluator;
  - the fraction of curator changes PASA-XG reproduced, which estimates recall;
  - changes not adopted, triaged with protein and peptide evidence (some may be real but
    not yet curated).
- This is the most direct evidence for "annotations that reflect current knowledge", and
  it doesn't depend on the truth set being complete.

### E5. GENCODE-XG: a population-scale protein-coding update (C5)

- **Candidates.** Download the Isopedia index and extract its full set of
  splice-junction chains (the README documents no export command, so this means reading
  the archive or using `isopedia-tools`). Isopedia mixes good signal with a lot of noise,
  so the filter needs its own study before we choose. Candidate criteria to explore:
  - recurrence across samples, tissues and studies;
  - support from both platforms (PacBio and ONT);
  - read counts and TPM per sample;
  - canonical splice motifs and junction support;
  - no evidence of internal priming or a truncated 5' end;
  - agreement with LRAA isoforms built from a subset of the source samples.
  Isoform ends aren't stored in the index (only junction chains), so the ends come from
  the read-end ranges Isopedia reports, or from our own LRAA runs.
- **Update:** PASA-XG on current GENCODE (v49) with defaults (protection, 40% coding,
  NMD rule).
- **Outputs:**
  - genes with a changed protein, new protein isoforms and new coding genes;
  - classes of protein change: N-terminal extension, alternative C-terminus, internal
    in-frame insertion or deletion, domain gain or loss;
  - expansion of distinct protein sequences vs GENCODE.
- **Orthogonal validation:**
  - peptides unique to new isoforms (large human proteomes, PeptideAtlas/PRIDE builds);
  - GENCODE Ribo-seq ORFs;
  - PhyloCSF/conservation;
  - Pfam integrity;
  - independent support in RefSeq and CHESS 3;
  - structure: ColabFold/ESMFold pLDDT, and Foldseek comparison against the canonical
    isoform, as in Sommer et al. (GPU needed for a few hundred examples).
- **Case studies:** a few loci with peptide support plus confident structure.

### E6 (optional). A non-model genome

Manatee: the LRGASP PacBio/ONT data used by Paniagua et al., applied as an update of
NCBI's annotation (Release 102). Scored with BUSCO and OMArk, and compared with their
long-read-driven annotation. This is a direct, cheap head-to-head with the most closely
related published strategy.

## 5. Engineering prerequisites

0. **Isopedia extraction:** dump the junction chains and per-sample counts from the
   downloaded index; build a quality-feature table for the filter study in E5.

1. **Chained updates.** A PASA-XG output GFF3 is the next run's `--annotation`. IDs
   already continue (`Namer` resumes `PXG_gene_<N>` / `.pxg<N>`); provenance must carry
   across runs (per-model attribute: input set, run, isoform IDs).
2. **Multiple samples in one run:** per-sample LRAA, then `merge_LRAA_GTFs.py` (with its
   provenance tracking file) inside the assemble stage.
3. **ESTs:** test the transcripts path with real ESTs (chaining, single-pass fragments).
   5'-partial ORFs take the EST path.
4. **Evaluation kit:** gffcompare wrappers for the views in section 3, compleasm/BUSCO,
   OMArk, protein diff vs the starting annotation, damage report, provenance summaries.
5. **Change report:** per-gene changelog (old vs new protein, class of change, evidence).
   It is also the release note for an annotation maintainer.
6. **Scale:** everything runs on the local host (28 cores, about 190 GB RAM, shared with
   jobs that mostly run niced). UGER is off-limits. The WDL stays secondary.

## 6. Compute plan

All local: 28 cores, about 190 GB RAM. LRAA runs per tissue and per sample are the main
cost, and they go in parallel within that budget.

- **Disk:** `/broad/hptmp` has a 3 TB quota shared with other projects (about 990 GB free on 2026-10-26). The Isopedia index needs about 412 GB
  (107 GB compressed plus 305 GB unpacked); delete the archive after unpacking. Keep raw
  read sets only as long as needed.
- **Structure predictions** need a GPU. They are limited to a few hundred case-study
  proteins, so a GPU host or ColabFold service only at the end; ESMFold on CPU is a
  fallback for small proteins.

## 7. Risks

- **Reference incompleteness** makes true discoveries look like false positives.
  Mitigation: precision against the union of catalogs, plus E4 and orthogonal evidence.
- **Comparator setup.** BRAKER3 long-read mode is experimental. Use the published
  protocol and report its settings; offer the authors' review if possible.
- **ONT accuracy:** report PacBio and ONT separately (LRGASP and Paniagua both found ONT
  weaker).
- **Isopedia** has no export: candidates come from our own LRAA runs, scored by Isopedia.
- **TransDecoder on short isoforms and ESTs:** partial ORFs go down the EST path by
  design. Report how many models come from each path.
- **Historic data availability** (dbEST, old GENCODE/TAIR releases): all public; confirm
  formats early.

## 8. Phasing

1. **Engineering** (section 5, items 1–4), plus E0 and E1 on Arabidopsis: the first
   whole-genome result and an end-to-end shakedown.
2. **E3 (Arabidopsis cross-generation) and E2 (Arabidopsis tissues):** the core
   cross-generation figure.
3. **Human E1/E2 on chromosomes locally,** then genome-wide on Terra; E4 time travel.
4. **E5 GENCODE-XG** and case studies; E6 if time allows.

Possible figures:

1. Concept: PASA-XG as the update layer across generations.
2. E1: draft → updated, vs comparators.
3. E3: generations alone and layered.
4. E2: saturation and order independence.
5. E4: curator adoption.
6. E5: proteome expansion with validation and structures.
7. Damage and run time (supplement or panel).

## 9. Decisions (2026-10-26)

1. Compute: all local on the 28-core host; no Terra for now.
2. Arabidopsis carries the first end-to-end demonstration.
3. Core paper: E1 (draft to current vs other tools), E3 (cross-generation), E4 (time
   travel), E5 (GENCODE-XG proteome). E0, E2 and E6 are supporting or supplementary.
4. Isopedia: download the index and extract all splice patterns. Selecting the best
   subset of a noisy catalog is its own study, to be designed before choosing.
5. Alternative-splicing analysis stays out of scope; it can be added later.
