"""Three-way diff on chunk_0's peaks:
  A = FIMO 5.5.5 chunk_0 output (this branch's run)
  B = FIMO 5.1.1 chunk_0 output (my controlled re-run)
  C = 2023 GT fimo.tsv filtered to chunk_0's 6,765 peaks

If A ⊆ C and B ⊆ C, and B ≈ what C contained on chunk_0's peaks, then the
2023 GT was made with FIMO 5.1.1-or-near and we understand all 3 datasets.
If B is much smaller than C∩chunk_0_peaks, then the GT used different flags
or a different MEME version we haven't pinned.
"""

import os, sys

A = "/fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/results/tcell-bm-validate/fimo_chunks/chunk_0/fimo.tsv"
B = "/fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/results/tcell-bm-validate/fimo511_chunk_0/fimo.tsv"
GT = "/fh/fast/setty_m/user/yhuang2/cell_cell_communication/tcell_depleted_bone_marrow/data/fimo.tsv"
CHUNK_FA = "/fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/results/tcell-bm-validate/fimo_chunks/chunk_0.fa"


def chunk0_peaks():
    out = set()
    with open(CHUNK_FA) as f:
        for line in f:
            if line.startswith(">"):
                out.add(line[1:].strip().split()[0])
    return out


def read_per_peak(path, only_peaks):
    """Return dict: peak -> set of (motif, alt, start, stop, strand, score, pval)"""
    per_peak = {}
    n_total = 0; n_in = 0
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            cols = line.rstrip("\n").split("\t")
            if cols[0] == "motif_id": continue
            if len(cols) < 8: continue
            n_total += 1
            seq = cols[2]
            if seq not in only_peaks: continue
            n_in += 1
            hit = (cols[0], cols[1], cols[3], cols[4], cols[5], cols[6], cols[7])
            per_peak.setdefault(seq, set()).add(hit)
    return per_peak, n_total, n_in


peaks0 = chunk0_peaks()
print(f"chunk_0 peaks: {len(peaks0):,}", flush=True)

print("\nreading 5.5.5 chunk_0 ...", flush=True)
A_pp, A_total, _ = read_per_peak(A, peaks0)
print(f"  {A_total:,} hits in file; {sum(len(v) for v in A_pp.values()):,} unique hits on {len(A_pp):,} peaks")

print("\nreading 5.1.1 chunk_0 ...", flush=True)
B_pp, B_total, _ = read_per_peak(B, peaks0)
print(f"  {B_total:,} hits in file; {sum(len(v) for v in B_pp.values()):,} unique hits on {len(B_pp):,} peaks")

print(f"\nreading GT (5.7 GB; filtering to chunk_0's {len(peaks0):,} peaks) ...", flush=True)
GT_pp, GT_total, GT_in = read_per_peak(GT, peaks0)
print(f"  {GT_total:,} rows in GT total; {GT_in:,} on chunk_0's peaks; {sum(len(v) for v in GT_pp.values()):,} unique on {len(GT_pp):,} peaks")

# Pairwise overlap on hit-tuples
def to_set(pp):
    out = set()
    for peak, hits in pp.items():
        for h in hits:
            out.add((peak,) + h)
    return out

A_set = to_set(A_pp); B_set = to_set(B_pp); GT_set = to_set(GT_pp)
print(f"\nA (5.5.5): {len(A_set):,}")
print(f"B (5.1.1): {len(B_set):,}")
print(f"GT (chunk_0 peaks): {len(GT_set):,}")
print()
print(f"A ∩ GT:  {len(A_set & GT_set):,}")
print(f"B ∩ GT:  {len(B_set & GT_set):,}")
print(f"A ∩ B:   {len(A_set & B_set):,}")
print(f"A \\ B (in 5.5.5 not 5.1.1, same chunk_0 input): {len(A_set - B_set):,}")
print(f"B \\ A (in 5.1.1 not 5.5.5, same chunk_0 input): {len(B_set - A_set):,}")
print(f"GT \\ B (in GT but my 5.1.1 didn't reproduce):    {len(GT_set - B_set):,}")
print(f"GT \\ A (in GT but my 5.5.5 didn't reproduce):    {len(GT_set - A_set):,}")

# Per-peak hit counts as a sanity check (sample 5 peaks)
print("\nSample 5 peaks, hit counts:")
import random
sample = random.sample(sorted(peaks0), 5)
for p in sample:
    a = len(A_pp.get(p, set()))
    b = len(B_pp.get(p, set()))
    g = len(GT_pp.get(p, set()))
    print(f"  {p:40s}  A(5.5.5)={a:>4}  B(5.1.1)={b:>4}  GT={g:>4}")
