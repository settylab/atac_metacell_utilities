#!/bin/bash
set -euo pipefail
export SCRATCH_DIR="${SCRATCH_DIR:-/loc/scratch/51961252}"
mkdir -p "$SCRATCH_DIR"
echo "[$(date)] streaming bucket diff on $(hostname); SCRATCH=$SCRATCH_DIR"
# Copy script into shared FS path; the script reads NEW + GT directly.
cp /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/exploration/fimo-version-investigation/validate_fimo_stream.py "$SCRATCH_DIR/validate_fimo_stream.py"
/home/yhuang2/micromamba/envs/seacells/bin/python "$SCRATCH_DIR/validate_fimo_stream.py" 2>&1
echo "[$(date)] done"
