"""Run settings and the per-input-type presets."""

from __future__ import annotations

import os

import shlex
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional

# minimap2 settings per --input_type. -ub everywhere: minimap2 picks each spliced
# read's transcribed strand from its splice motifs (the ts tag), and LRAA resolves
# strands from ts on its own, so stranded and unstranded libraries need no option.
MINIMAP2_PRESETS: Dict[str, List[str]] = {
    "pacbio_hifi": ["-ax", "splice:hq", "-ub"],
    "ont": ["-ax", "splice", "-ub"],
    # PASA's settings for transcript alignments (process_minimap2_alignments.pl)
    "transcripts": ["-ax", "splice", "--secondary=no", "-O6,24", "-B4", "-L", "-cs", "-ub"],
}

LRAA_PRESETS: Dict[str, List[str]] = {
    "pacbio_hifi": ["--HiFi"],
    "ont": [],
    "transcripts": ["--pasa_mode", "--preassembled_input"],
}

INPUT_TYPES = tuple(MINIMAP2_PRESETS)


@dataclass
class Settings:
    genome: str
    input_type: str
    out_dir: str
    reads: List[str] = field(default_factory=list)
    bam: Optional[str] = None
    annotation: Optional[str] = None
    cpu: int = 4
    chain_fragments: bool = False
    minimap2_opts: Optional[str] = None
    max_intron_length: int = 100000
    lraa_cmd: str = "LRAA"
    lraa_opts: Optional[str] = None
    transdecoder_dir: Optional[str] = None
    genetic_code: str = "Universal"
    min_orf_aa: int = 100
    pfam_hits: Optional[str] = None
    blastp_hits: Optional[str] = None
    # annotation update
    max_rounds: int = 3
    run_label: str = ""                # recorded in pasa_xg_history (default: out dir name)
    protect_supported_models: bool = True
    additive: bool = False             # never replace, merge or split existing models
    trust_fl_status: bool = True       # PASA's --TRUST_FL_STATUS
    reject_nmd_targets: bool = True
    nmd_stop_to_junction: int = 50
    min_percent_coding: float = 40.0   # PASA's MIN_PERCENT_PROT_CODING
    min_fl_orf_aa: int = 0             # PASA's MIN_FL_ORF_SIZE (0: off)
    min_tpm: float = 0.0
    min_isoform_fraction: float = 0.0

    def __post_init__(self) -> None:
        if not self.run_label:
            self.run_label = os.path.basename(os.path.normpath(self.out_dir))
        if self.input_type not in INPUT_TYPES:
            raise ValueError(f"--input_type must be one of {', '.join(INPUT_TYPES)}")
        if bool(self.reads) == bool(self.bam):
            raise ValueError("give either reads/transcripts (--reads) or an aligned --bam, not both")

    def minimap2_args(self) -> List[str]:
        if self.minimap2_opts:
            return shlex.split(self.minimap2_opts)
        return list(MINIMAP2_PRESETS[self.input_type])

    def lraa_args(self) -> List[str]:
        args = list(LRAA_PRESETS[self.input_type])
        if self.chain_fragments and "--pasa_mode" not in args:
            args.append("--pasa_mode")
        if self.lraa_opts:
            args += shlex.split(self.lraa_opts)
        return args

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
