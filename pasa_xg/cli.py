"""pasa-xg command line."""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys

from . import __version__
from .config import INPUT_TYPES, Settings
from .pipeline import STAGES, run_all


def _run_parser(p: argparse.ArgumentParser) -> None:
    req = p.add_argument_group("inputs")
    req.add_argument("--genome", required=True, help="genome FASTA")
    src = req.add_mutually_exclusive_group(required=True)
    src.add_argument("--reads", nargs="+", default=[],
                     help="reads or transcripts, FASTA or FASTQ (gzip ok), aligned with minimap2")
    src.add_argument("--bam", help="coordinate-sorted alignments instead of --reads (skips alignment)")
    src.add_argument("--isoforms_gtf", help="finished genome-aligned isoforms (e.g. a merged catalog), used "
                                            "as they are: skips alignment and assembly, so nothing is "
                                            "chained. TSS and PolyA attributes are kept")
    req.add_argument("--input_type", choices=INPUT_TYPES,
                     help="pacbio_hifi, ont, or transcripts (Trinity/StringTie assemblies, ESTs, "
                          "cDNAs). Sets minimap2's and LRAA's options together. Required with "
                          "--reads or --bam.")
    req.add_argument("--annotation", help="gene annotation (GTF or GFF3). Guides LRAA; with CDS "
                                          "features it is also compared and updated.")
    req.add_argument("--out", required=True, dest="out_dir", help="run directory")
    p.add_argument("--cpu", type=int, default=4)
    p.add_argument("--from", dest="force_from", choices=STAGES,
                   help="rerun this stage and everything after it")
    p.add_argument("--stop_after", choices=STAGES, help="stop after this stage")

    a = p.add_argument_group("alignment and assembly")
    a.add_argument("--chain_fragments", action="store_true",
                   help="reads only: chain fragmented reads that share introns (LRAA --pasa_mode); "
                        "always on for --input_type transcripts")
    a.add_argument("--minimap2_opts", help="replace the preset minimap2 options (quoted string)")
    a.add_argument("--max_intron_length", type=int, default=100000)
    a.add_argument("--lraa_cmd", default=os.environ.get("PASA_XG_LRAA", "LRAA"),
                   help="LRAA command (default: $PASA_XG_LRAA or LRAA)")
    a.add_argument("--lraa_opts", help="extra LRAA options (quoted string)")

    c = p.add_argument_group("coding regions")
    c.add_argument("--transdecoder_cmd", default="pytransdecoder",
                   help="pytransdecoder command (PyPI pytransdecoder >= 0.2.1)")
    c.add_argument("--genetic_code", default="Universal")
    c.add_argument("--min_orf_aa", type=int, default=100)
    c.add_argument("--pfam_hits", help="hmmscan --domtblout of longest_orfs.pep against Pfam")
    c.add_argument("--blastp_hits", help="blastp -outfmt 6 of longest_orfs.pep")

    u = p.add_argument_group("annotation update")
    u.add_argument("--max_rounds", type=int, default=3)
    u.add_argument("--run_label", default="",
                   help="name of this update, recorded in each changed model's "
                        "pasa_xg_history (default: the --out directory name). Feed this "
                        "run's annotation.updated.gff3 to a later run as --annotation to "
                        "keep updating it; histories accumulate.")
    u.add_argument("--allow_nmd_targets", dest="reject_nmd_targets", action="store_false",
                   help="let updates produce models that look like NMD targets. By default "
                        "a model whose stop codon lies more than --nmd_stop_to_junction "
                        "bases upstream of its last splice junction is rejected.")
    u.add_argument("--nmd_stop_to_junction", type=int, default=50)
    u.add_argument("--min_percent_coding", type=float, default=40.0,
                   help="full-length isoforms must be at least this %% coding (PASA's "
                        "MIN_PERCENT_PROT_CODING)")
    u.add_argument("--min_fl_orf_aa", type=int, default=0,
                   help="a full-length isoform whose protein is longer than this passes "
                        "regardless of --min_percent_coding (PASA's MIN_FL_ORF_SIZE; 0: off)")
    u.add_argument("--no_trust_fl_status", dest="trust_fl_status", action="store_false",
                   help="PASA's default: a compatible full-length isoform that the gene model "
                        "extends past is treated as an EST. By default (PASA's "
                        "--TRUST_FL_STATUS) it stays full-length, since LRAA trims isoform ends.")
    u.add_argument("--additive", action="store_true",
                   help="for curated annotations: never replace, restitch or split an "
                        "existing model, even one the reads don't support (absence of "
                        "evidence in these samples isn't evidence against it). Conflicting "
                        "structures are added as alternative isoforms or new genes; UTR "
                        "extensions, compatible CDS refinements and guarded gene merges "
                        "still apply. Default: PASA's behavior, which replaces unsupported "
                        "models (suits drafts).")
    u.add_argument("--replace_supported_models", dest="protect_supported_models",
                   action="store_false",
                   help="PASA's behavior: let an incompatible isoform replace (or be stitched "
                        "into) a model even when other isoforms support that model's own "
                        "structure. By default such a model is kept and the new structure "
                        "is added as an alternative isoform.")
    u.add_argument("--no_reference_frame_check", dest="reference_frame_check", action="store_false",
                   help="keep every isoform's ORF as called. By default, an isoform holding the gene's start "
                        "codon upstream of its ORF is read from that start (noncoding when that ORF is short or "
                        "an NMD target), and an isoform whose 5' or 3' end falls inside the reference CDS "
                        "without TSS / polyA support is read in the reference frame as a partial "
                        "(cds_start_NF / cds_end_NF) instead of from a downstream Met or a new stop")
    u.add_argument("--min_sample_fraction", type=float, default=0.01,
                   help="noncoding pass, for isoforms carrying n_samples and n_samples_total (multi-sample "
                        "merged catalogs): the evidence filter's minimum of 2 samples becomes "
                        "max(2, fraction x n_samples_total), and its other sample counts move up with it")
    u.add_argument("--coding_only", dest="include_noncoding", action="store_false",
                   help="skip the noncoding pass and update only coding isoforms, as PASA does. By "
                        "default, after the coding rounds, isoforms the coding update turned away "
                        "that have a new structure are added: to their gene, labeled pasa_xg_coding "
                        "low / nmd / none, or as novel noncoding genes, and every gene is labeled "
                        "pasa_xg_gene_class coding / likely_noncoding. Isoforms carrying an "
                        "n_samples GTF attribute (multi-sample merged catalogs) must also pass an "
                        "evidence filter (docs/NONCODING_DESIGN.md)")
    u.add_argument("--min_tpm", type=float, default=0.0,
                   help="isoforms below this TPM do not update the annotation")
    u.add_argument("--min_isoform_fraction", type=float, default=0.0,
                   help="isoforms below this fraction of their gene do not update the annotation")

    t = p.add_argument_group("TSS and polyA sites (docs/TSS_POLYA_DESIGN.md)")
    t.add_argument("--tss_sites", nargs="+", default=[],
                   help="known TSS sites, BED with strand (e.g. FANTOM5 CAGE peaks, refTSS). "
                        "With these, or LRAA's own TSS/PolyA calls on the isoforms, an isoform "
                        "end only moves a model end that no site supports, and models are "
                        "never shortened")
    t.add_argument("--polya_sites", nargs="+", default=[],
                   help="known polyA sites, BED with strand (e.g. PolyASite clusters)")
    t.add_argument("--site_window", type=int, default=50,
                   help="an end within this many bases of a site is at that site (LRAA's "
                        "max_dist_between_alt_TSS_sites)")
    t.add_argument("--alt_end_isoforms", choices=("protein", "all"), default="protein",
                   help="protein: an isoform that differs from a model only at its ends is "
                        "added only when its protein differs (and, with sites, only when those "
                        "ends are supported); its supported ends are recorded on the model as "
                        "alt_tss_sites / alt_polya_sites. all: PASA's behavior")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="pasa-xg", description="PASA-XG: transcript-evidence "
                                "genome annotation across sequencing generations")
    p.add_argument("--version", action="version", version=f"PASA-XG {__version__}")
    sub = p.add_subparsers(dest="command", required=True)
    _run_parser(sub.add_parser("run", help="run the pipeline (resumes a stopped run)"))
    st = sub.add_parser("status", help="stage states of a run")
    st.add_argument("--out", required=True, dest="out_dir")
    args = p.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")
    if args.command == "status":
        with open(os.path.join(args.out_dir, "manifest.json")) as fh:
            m = json.load(fh)
        for name in STAGES:
            st = m["stages"].get(name)
            print(f"{name:14s} {st['status'] if st else 'not run'}"
                  + (f"  ({st.get('seconds', '?')} s)" if st and "seconds" in st else ""))
        return 0

    opts = {k: v for k, v in vars(args).items()
            if k not in ("command", "force_from", "stop_after")}
    try:
        s = Settings(**opts)
    except ValueError as e:
        p.error(str(e))
    out = run_all(s, force_from=args.force_from, stop_after=args.stop_after)
    summary = {stage: {k: v for k, v in o.items() if not k.endswith("_path")} for stage, o in out.items()}
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
