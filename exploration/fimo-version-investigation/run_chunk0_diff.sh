#!/bin/bash
set -euo pipefail
echo "[$(date)] three-way chunk_0 diff on $(hostname)"
/home/yhuang2/micromamba/envs/seacells/bin/python \
  /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/exploration/fimo-version-investigation/diff_chunk0_three_way.py 2>&1
echo "[$(date)] done"
