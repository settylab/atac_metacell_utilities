#!/bin/bash
# Hypothesis: GT was made with FIMO 5.1.1 DEFAULT flags (no --no-qvalue);
# my earlier 5.1.1 probe used --no-qvalue which may filter differently in
# FIMO 5.1.1 (vs 5.5.5 where --no-qvalue only skips q-value computation).
# Re-run 5.1.1 on chunk_0 with NO flags except --oc, matching the GT
# command-line `fimo -oc <dir> <meme> <fasta>` exactly.
set -euo pipefail
cd /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s
FIMO=/app/software/MEME/5.1.1-foss-2019b-Perl-5.30.0-Python-3.7.4/bin/fimo

OUT=results/tcell-bm-validate/fimo511_noflags_chunk_0
LOG=results/tcell-bm-validate/logs/fimo511_noflags_chunk_0.out
mkdir -p "$OUT" "$(dirname "$LOG")"

START=$(date +%s)
echo "[$(date)] fimo 5.1.1 NO-FLAGS probe on $(hostname): $($FIMO --version)"
"$FIMO" -oc "$OUT" \
     data/cis-bp-tf-information.meme \
     results/tcell-bm-validate/fimo_chunks/chunk_0.fa \
     2> "$LOG"
rm -f "$OUT/cisml.xml"
END=$(date +%s)
echo "[$(date)] done; elapsed $((END - START))s"
echo "fimo.tsv hit count: $(($(wc -l < "$OUT/fimo.tsv") - 1))"
