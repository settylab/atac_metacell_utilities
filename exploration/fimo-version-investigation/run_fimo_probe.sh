#!/bin/bash
set -euo pipefail
cd /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s

CHUNK_DIR="results/tcell-bm-validate/fimo_chunks/chunk_0"
LOG="results/tcell-bm-validate/logs/fimo_chunk_0.out"
mkdir -p "$CHUNK_DIR" "$(dirname "$LOG")"

# Use seacells env's fimo (5.5.5) — already on PATH for srun-spawned shell
# only if we prepend; pass the absolute path instead.
FIMO=/home/yhuang2/micromamba/envs/seacells/bin/fimo

START=$(date +%s)
echo "[$(date)] fimo probe on $(hostname); chunk_0 = $(grep -c '^>' results/tcell-bm-validate/fimo_chunks/chunk_0.fa) sequences"

# /usr/bin/time -v gives us peak RSS for resource sizing.
/usr/bin/time -v "$FIMO" --no-pgc --thresh 1e-4 --no-qvalue \
    -oc "$CHUNK_DIR" \
    data/cis-bp-tf-information.meme \
    results/tcell-bm-validate/fimo_chunks/chunk_0.fa \
    2> "$LOG"
rm -f "$CHUNK_DIR/cisml.xml"

END=$(date +%s)
echo "[$(date)] fimo probe done; elapsed $((END - START))s"
echo "fimo.tsv hit count: $(($(wc -l < "$CHUNK_DIR/fimo.tsv") - 1))"
echo
echo "=== /usr/bin/time -v summary (last 25 lines of $LOG) ==="
tail -25 "$LOG"
