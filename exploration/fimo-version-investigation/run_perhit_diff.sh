#!/bin/bash
set -euo pipefail
export SCRATCH_DIR="${SCRATCH_DIR:-/loc/scratch/51961252}"
mkdir -p "$SCRATCH_DIR"
echo "[$(date)] per-hit bucket diff on $(hostname); SCRATCH=$SCRATCH_DIR"
/home/yhuang2/micromamba/envs/seacells/bin/python \
  /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/exploration/fimo-version-investigation/validate_fimo_perhit.py 2>&1
echo "[$(date)] done"
