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
    req.add_argument("--input_type", required=True, choices=INPUT_TYPES,
                     help="pacbio_hifi, ont, or transcripts (Trinity/StringTie assemblies, ESTs, "
                          "cDNAs). Sets minimap2's and LRAA's options together.")
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
    c.add_argument("--transdecoder_dir", help="directory holding TransDecoder.LongOrfs/Predict")
    c.add_argument("--genetic_code", default="Universal")
    c.add_argument("--min_orf_aa", type=int, default=100)
    c.add_argument("--pfam_hits", help="hmmscan --domtblout of longest_orfs.pep against Pfam")
    c.add_argument("--blastp_hits", help="blastp -outfmt 6 of longest_orfs.pep")

    u = p.add_argument_group("annotation update")
    u.add_argument("--max_rounds", type=int, default=3)
    u.add_argument("--replace_supported_models", dest="protect_supported_models",
                   action="store_false",
                   help="PASA's behavior: let an incompatible isoform replace (or be stitched "
                        "into) a model even when other isoforms support that model's own "
                        "structure. By default such a model is kept and the new structure "
                        "is added as an alternative isoform.")
    u.add_argument("--min_tpm", type=float, default=0.0,
                   help="isoforms below this TPM do not update the annotation")
    u.add_argument("--min_isoform_fraction", type=float, default=0.0,
                   help="isoforms below this fraction of their gene do not update the annotation")


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
