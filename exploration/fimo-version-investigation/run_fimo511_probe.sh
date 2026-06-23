#!/bin/bash
# Hypothesis test: re-run chunk_0 with FIMO 5.1.1 (the version that produced
# the 2023 ground truth) using identical inputs and flags as the 5.5.5 run.
# If 5.1.1 output exactly matches the 2023 GT fimo.tsv on chunk_0's peaks,
# then the discrepancy vs our 5.5.5 run is purely FIMO version drift.
set -euo pipefail
cd /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s

FIMO=/app/software/MEME/5.1.1-foss-2019b-Perl-5.30.0-Python-3.7.4/bin/fimo

echo "[$(date)] fimo 5.1.1 probe on $(hostname): $($FIMO --version)"
echo "input: $(grep -c '^>' results/tcell-bm-validate/fimo_chunks/chunk_0.fa) sequences"

OUT=results/tcell-bm-validate/fimo511_chunk_0
LOG=results/tcell-bm-validate/logs/fimo511_chunk_0.out
mkdir -p "$OUT" "$(dirname "$LOG")"

START=$(date +%s)
"$FIMO" --thresh 1e-4 --no-qvalue \
     --oc "$OUT" \
     data/cis-bp-tf-information.meme \
     results/tcell-bm-validate/fimo_chunks/chunk_0.fa \
     2> "$LOG"
rm -f "$OUT/cisml.xml"
END=$(date +%s)
echo "[$(date)] done; elapsed $((END - START))s"
echo "fimo.tsv hit count: $(($(wc -l < "$OUT/fimo.tsv") - 1))"
echo
echo "=== FIMO 5.1.1 first 3 rows ==="
head -4 "$OUT/fimo.tsv"
echo "=== FIMO 5.5.5 first 3 rows for comparison ==="
head -4 results/tcell-bm-validate/fimo_chunks/chunk_0/fimo.tsv
