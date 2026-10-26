#!/bin/bash
# Development LRAA: the pinned pasa_mode_v2 checkout, run in the LRAA core image
# (python deps), with the project tree bound. Use as --lraa_cmd or $PASA_XG_LRAA.
P=/broad/hptmp/bhaas/PASAPIPELINE_PAPER
SIF=${LRAA_SIF:-/broad/hptmp/bhaas/resources/containers/lraa-core_0.44.3-5c002da.sif}
mkdir -p $P/PASA-XG/scratch/tmp
SRC=$P/LRAA_pasa_mode/LongReadAlignmentAssembler
if [ "$1" == "--version" ]; then  # include the commit: PASA-XG keys its cache on this
  echo "LRAA $(grep -m1 -oP '^VERSION = "\K[^"]+' $SRC/LRAA) $(git -C $SRC rev-parse --short HEAD) $(git -C $SRC rev-parse --abbrev-ref HEAD)"
  exit 0
fi
exec /broad/hptmp/bhaas/tools/apptainer-1.5.3/bin/apptainer exec --containall --no-mount hostfs \
  --pwd "$PWD" -B $P --env TMPDIR=$P/PASA-XG/scratch/tmp,OMP_NUM_THREADS=1 $SIF \
  $SRC/LRAA "$@"
