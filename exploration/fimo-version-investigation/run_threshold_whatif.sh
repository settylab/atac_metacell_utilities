#!/bin/bash
set -euo pipefail
echo "[$(date)] threshold what-if on $(hostname)"
/home/yhuang2/micromamba/envs/seacells/bin/python \
  /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/exploration/fimo-version-investigation/threshold_what_if.py 2>&1
echo "[$(date)] done"
