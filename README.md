# PASA-XG

Transcript-evidence genome annotation across sequencing generations: ESTs and
cDNAs, short-read assemblies (Trinity, StringTie) and PacBio/ONT long reads.
A rebuild of the [PASA](https://github.com/PASApipeline/PASApipeline) pipeline:
minimap2 aligns, [LRAA](https://github.com/MethodsDev/LongReadAlignmentAssembler)
assembles and quantifies isoforms, TransDecoder finds coding regions, and PASA's
annotation comparison and update logic runs on the result.

Status: early development. See [DESIGN.md](DESIGN.md) and [STATUS.md](STATUS.md).

## Quick start

```bash
pip install -e .            # plus minimap2, samtools, LRAA, TransDecoder on PATH
pasa-xg run --genome genome.fa --reads reads.fq.gz --input_type pacbio_hifi \
    --annotation annot.gff3 --cpu 16 --out run/
pasa-xg status --out run/
```

`--input_type` is `pacbio_hifi`, `ont` or `transcripts`; `--bam` replaces `--reads`
to start from alignments. A stopped run resumes where it left off; `--from <stage>`
reruns a stage and everything after it.

## Development

- LRAA: branch `pasa_mode_v2` (commit 2a6524d) of MethodsDev/LongReadAlignmentAssembler,
  which adds `--pasa_mode` and `--preassembled_input`. `scripts/lraa_dev.sh` runs that
  checkout in the LRAA image; point `$PASA_XG_LRAA` at it.
- Tests: `pytest tests`.
