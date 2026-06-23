"""What-if analysis on tightening --thresh for FIMO 5.5.5.

For chunk_0 (6,765 peaks), distribute hits by p-value across thresholds
{1e-4, 5e-5, 1e-5, 5e-6, 1e-6}. Separately for:
  - NEW chunk_0 hits that match a GT hit at the same position (these are
    the hits we MUST preserve to keep A ⊇ GT)
  - NEW chunk_0 hits that are NEW-only (extras at the threshold edge)

Output: a table showing for each candidate threshold, how many shared
(GT-matching) hits survive and how many NEW-only hits survive.
"""

import sys

NEW_chunk0 = "/fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/results/tcell-bm-validate/fimo_chunks/chunk_0/fimo.tsv"
GT_full   = "/fh/fast/setty_m/user/yhuang2/cell_cell_communication/tcell_depleted_bone_marrow/data/fimo.tsv"
CHUNK_FA  = "/fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/results/tcell-bm-validate/fimo_chunks/chunk_0.fa"


def chunk0_peaks():
    out = set()
    with open(CHUNK_FA) as f:
        for line in f:
            if line.startswith(">"):
                out.add(line[1:].strip().split()[0])
    return out


def read_hits(path, peak_filter=None):
    """Yield (motif, alt, seq, start, stop, strand) -> p-value (float).
    Returns dict; if peak_filter is given, only includes hits on those peaks."""
    hits = {}
    n_total = 0
    n_kept = 0
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[0] == "motif_id":
                continue
            if len(cols) < 8:
                continue
            n_total += 1
            seq = cols[2]
            if peak_filter is not None and seq not in peak_filter:
                continue
            n_kept += 1
            try:
                pv = float(cols[7])
            except ValueError:
                continue
            key = (cols[0], cols[1], cols[2], cols[3], cols[4], cols[5])
            hits[key] = pv
    print(f"  {path}: {n_total:,} rows total, {n_kept:,} on chunk_0 peaks, {len(hits):,} unique hits", flush=True)
    return hits


peaks0 = chunk0_peaks()
print(f"chunk_0 peaks: {len(peaks0):,}", flush=True)

print("\nreading NEW chunk_0 ...", flush=True)
new_hits = read_hits(NEW_chunk0)

print("\nreading GT (filtered to chunk_0 peaks) ...", flush=True)
gt_hits = read_hits(GT_full, peak_filter=peaks0)

# Partition NEW into shared-with-GT and NEW-only
shared_keys = set(new_hits) & set(gt_hits)
new_only_keys = set(new_hits) - set(gt_hits)
print(f"\nNEW chunk_0 total: {len(new_hits):,}")
print(f"  shared with GT: {len(shared_keys):,}")
print(f"  NEW-only:       {len(new_only_keys):,}")

# Threshold sweep
thresholds = [1e-4, 5e-5, 1e-5, 5e-6, 1e-6]
print()
print(f"  {'threshold':>10s}  {'shared survive':>14s}  {'NEW-only survive':>17s}  {'GT preserved?':>14s}")
print(f"  {'-'*10:>10s}  {'-'*14:>14s}  {'-'*17:>17s}  {'-'*14:>14s}")
shared_pvs = [new_hits[k] for k in shared_keys]
new_only_pvs = [new_hits[k] for k in new_only_keys]
N_shared = len(shared_pvs)
N_new_only = len(new_only_pvs)
for t in thresholds:
    n_shared_surv = sum(1 for p in shared_pvs if p < t)
    n_new_only_surv = sum(1 for p in new_only_pvs if p < t)
    n_shared_lost = N_shared - n_shared_surv
    preserved = "YES" if n_shared_lost == 0 else f"NO (lose {n_shared_lost:,})"
    print(f"  {t:10.0e}  {n_shared_surv:>14,}  {n_new_only_surv:>17,}  {preserved:>14s}")

# Distribution of p-values among shared (GT) hits to understand the curve
print()
print("Distribution of p-values among shared (GT-matching) hits, NEW chunk_0:")
import math
edges = [(-30,-9), (-9,-7), (-7,-6), (-6,-5), (-5,-4.5), (-4.5,-4.2), (-4.2,-4.1), (-4.1,-4.0)]
for lo, hi in edges:
    n = sum(1 for p in shared_pvs if p > 0 and lo <= math.log10(p) < hi)
    bar = "#" * max(1, int(40 * n / max(N_shared, 1))) if n else ""
    print(f"  p in [10^{lo:>+5.1f}, 10^{hi:>+4.1f}):  {n:>10,}  {bar}")
