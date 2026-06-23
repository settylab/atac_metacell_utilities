#!/bin/bash
set -euo pipefail
cd /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s
START=$(date +%s)
echo "[$(date)] seq_gl R run on $(hostname)"
/home/yhuang2/micromamba/envs/r_env/bin/Rscript \
  results/tcell-bm-validate/seq_gl_for_diff.R \
  results/tcell-bm-validate/peaks.bed \
  results/tcell-bm-validate/all_seqs.r.fa \
  150 hg38 2>&1
END=$(date +%s)
echo "[$(date)] done; elapsed $((END - START))s"
ls -la results/tcell-bm-validate/all_seqs.r.fa
