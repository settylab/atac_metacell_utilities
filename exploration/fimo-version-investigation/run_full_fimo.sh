#!/bin/bash
# Full N=32 FIMO scatter via snakemake inside an srun on gizmok124.
# srun --ntasks=1 to prevent the parent allocation's --ntasks default
# from launching multiple copies of this script.
set -euo pipefail
cd /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s
export PATH=/home/yhuang2/micromamba/envs/seacells/bin:$PATH
export PYTHONPATH=/tmp/pyfaidx_validate
START=$(date +%s)
echo "[$(date)] full FIMO scatter on $(hostname); 32 chunks via snakemake --cores 32 --until fimo"
echo "[$(date)] FIMO version: $(fimo 2>&1 | head -1)"
/home/yhuang2/micromamba/envs/omnilib/bin/snakemake \
  --configfile results/tcell-bm-validate/config.yaml \
  --cores 32 \
  --until fimo \
  --rerun-incomplete \
  --keep-going 2>&1 | tee results/tcell-bm-validate/snakemake_full.log
END=$(date +%s)
echo "[$(date)] full run done; elapsed $((END - START))s"
ls -la results/tcell-bm-validate/fimo_result/fimo.tsv
echo "fimo.tsv hit count: $(($(wc -l < results/tcell-bm-validate/fimo_result/fimo.tsv) - 1))"
