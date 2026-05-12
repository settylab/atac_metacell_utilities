#!/bin/bash
# Test whether scatter introduces semantic drift via per-chunk vs global background.
# 1. Compute 0th-order Markov bg from FULL all_seqs.fa.
# 2. Re-run FIMO 5.5.5 on chunk_0.fa with --bfile pointing at the global bg.
# 3. Hit count + sample compared against the existing scatter chunk_0 output
#    (which used FIMO's default per-chunk bg).
set -euo pipefail
cd /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s

echo "[$(date)] bg-test on $(hostname)"
which fasta-get-markov || ls /home/yhuang2/micromamba/envs/seacells/bin/fasta-get-markov 2>&1
FIMO=/home/yhuang2/micromamba/envs/seacells/bin/fimo
FGM=/home/yhuang2/micromamba/envs/seacells/bin/fasta-get-markov

BG=results/tcell-bm-validate/all_seqs.bg
OUT=results/tcell-bm-validate/fimo_chunk0_globalbg
LOG=results/tcell-bm-validate/logs/fimo_chunk0_globalbg.out
mkdir -p "$OUT" "$(dirname "$LOG")"

if [ ! -s "$BG" ]; then
  echo "[$(date)] computing 0th-order Markov bg from full all_seqs.fa (216,477 records, 32.5 MB)"
  "$FGM" -m 0 results/tcell-bm-validate/all_seqs.fa "$BG"
fi
echo "=== global bg (from all_seqs.fa) ==="; cat "$BG"

START=$(date +%s)
echo "[$(date)] re-running FIMO 5.5.5 on chunk_0 with --bfile <global-bg>"
"$FIMO" --bfile "$BG" --no-pgc --thresh 1e-4 --no-qvalue \
    --oc "$OUT" \
    data/cis-bp-tf-information.meme \
    results/tcell-bm-validate/fimo_chunks/chunk_0.fa \
    2> "$LOG"
rm -f "$OUT/cisml.xml" "$OUT/fimo.html"
END=$(date +%s)
echo "[$(date)] done; elapsed $((END - START))s"
echo "hits with global bg: $(($(wc -l < "$OUT/fimo.tsv") - 1))"
echo "hits with per-chunk bg (existing scatter run): $(($(wc -l < results/tcell-bm-validate/fimo_chunks/chunk_0/fimo.tsv) - 1))"
