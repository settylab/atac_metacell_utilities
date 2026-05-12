"""Per-hit regression diff: NEW vs ground-truth FIMO TSVs, keyed by the
FULL position tuple (motif_id, motif_alt_id, sequence_name, start, stop,
strand) so multi-hit peaks aren't collapsed.

Streams both files in a single pass into hash-bucketed temp files keyed
by sequence_name, then compares bucket-by-bucket to keep RSS low.

Output (per-bucket and aggregated over all buckets):
- shared / NEW-only / GT-only hit counts (true per-hit, not per-triple)
- on shared hits: score match (byte-equal string), score match (np.isclose)
- on shared hits: p-value match (byte-equal), p-value match (np.isclose)
- distribution (log10 |Delta|) of p-value drift on shared hits
- distribution of NEW-only hits: are they at the p-value/score boundary?
"""

import os, sys, hashlib, tempfile, collections, math

NEW = "/fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/results/tcell-bm-validate/fimo_result/fimo.tsv"
GT  = "/fh/fast/setty_m/user/yhuang2/cell_cell_communication/tcell_depleted_bone_marrow/data/fimo.tsv"
N_BUCKETS = 16


def bucket_of(name):
    return int(hashlib.md5(name.encode()).hexdigest()[:4], 16) % N_BUCKETS


def split_to_buckets(path, prefix):
    """Stream rows into N_BUCKETS files keyed by sequence_name hash.
    Output cols: motif_id, motif_alt_id, sequence_name, start, stop, strand, score, pvalue
    (we drop matched_sequence column to shrink temp files)
    """
    handles = [open(f"{prefix}.{b}.tsv", "w") for b in range(N_BUCKETS)]
    n = 0
    try:
        with open(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if cols[0] == "motif_id":
                    continue
                if len(cols) < 8:
                    continue
                b = bucket_of(cols[2])
                # motif_id  motif_alt_id  seq  start  stop  strand  score  pval
                handles[b].write(
                    f"{cols[0]}\t{cols[1]}\t{cols[2]}\t{cols[3]}\t{cols[4]}\t{cols[5]}\t{cols[6]}\t{cols[7]}\n"
                )
                n += 1
    finally:
        for h in handles:
            h.close()
    return n


def _safe_log10_abs(x):
    if x <= 0:
        return None
    return math.log10(x)


def compare_bucket(new_path, gt_path):
    """Per-hit comparison for one bucket."""
    def read(p):
        hits = {}    # (motif, alt, seq, start, stop, strand) -> (score_str, pval_str)
        with open(p) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 8:
                    continue
                key = (cols[0], cols[1], cols[2], cols[3], cols[4], cols[5])
                # If duplicate keys (FIMO shouldn't emit them, but be safe), keep first
                if key not in hits:
                    hits[key] = (cols[6], cols[7])
        return hits

    new_h = read(new_path)
    gt_h = read(gt_path)

    new_keys = set(new_h)
    gt_keys = set(gt_h)
    shared = new_keys & gt_keys
    new_only = new_keys - gt_keys
    gt_only = gt_keys - new_keys

    # On shared hits: score and p-value equality
    score_eq_str = 0
    score_eq_close = 0
    pval_eq_str = 0
    pval_eq_close = 0
    pval_log10_diffs = []  # log10 |new_p - gt_p| for non-equal
    score_diffs = []
    for k in shared:
        new_score, new_pv = new_h[k]
        gt_score, gt_pv = gt_h[k]
        if new_score == gt_score:
            score_eq_str += 1
        try:
            ns, gs = float(new_score), float(gt_score)
            if abs(ns - gs) <= 1e-6 * max(abs(ns), abs(gs), 1.0):
                score_eq_close += 1
            else:
                score_diffs.append(abs(ns - gs))
        except ValueError:
            pass
        if new_pv == gt_pv:
            pval_eq_str += 1
        try:
            npv, gpv = float(new_pv), float(gt_pv)
            if abs(npv - gpv) <= 1e-6 * max(abs(npv), abs(gpv), 1.0):
                pval_eq_close += 1
            else:
                d = abs(npv - gpv)
                lg = _safe_log10_abs(d)
                if lg is not None:
                    pval_log10_diffs.append(lg)
        except ValueError:
            pass

    # NEW-only: where do they sit on p-value scale? Threshold was 1e-4 (= -4 in log10).
    new_only_pv_log10 = []
    new_only_score = []
    for k in new_only:
        s, p = new_h[k]
        try:
            new_only_score.append(float(s))
            lg = _safe_log10_abs(float(p))
            if lg is not None:
                new_only_pv_log10.append(lg)
        except ValueError:
            pass

    return {
        "new_hits": len(new_keys),
        "gt_hits": len(gt_keys),
        "shared": len(shared),
        "new_only": len(new_only),
        "gt_only": len(gt_only),
        "score_eq_str": score_eq_str,
        "score_eq_close": score_eq_close,
        "pval_eq_str": pval_eq_str,
        "pval_eq_close": pval_eq_close,
        "pval_log10_diffs": pval_log10_diffs,
        "score_diffs": score_diffs,
        "new_only_pv_log10": new_only_pv_log10,
        "new_only_score": new_only_score,
    }


def histogram(values, bin_edges, title):
    """Print an ASCII histogram (bin label, count, bar) over the supplied bin_edges."""
    counts = [0] * (len(bin_edges) - 1)
    for v in values:
        for i in range(len(bin_edges) - 1):
            if bin_edges[i] <= v < bin_edges[i + 1]:
                counts[i] += 1
                break
        else:
            counts[-1] += 1
    mx = max(counts) if counts else 1
    print(f"\n{title}  (n={len(values):,})")
    for i in range(len(counts)):
        lo, hi = bin_edges[i], bin_edges[i + 1]
        bar = "#" * max(1, int(40 * counts[i] / mx)) if counts[i] else ""
        print(f"  [{lo:+6.1f}, {hi:+6.1f})  {counts[i]:>12,}  {bar}")


def main():
    scratch = os.environ.get("SCRATCH_DIR", "/tmp")
    with tempfile.TemporaryDirectory(prefix="fimo_diff_perhit_", dir=scratch) as td:
        print(f"bucketing into {td}/", flush=True)
        print(f"  pass 1: NEW", flush=True)
        n_new = split_to_buckets(NEW, f"{td}/new")
        print(f"    {n_new:,} rows", flush=True)
        print(f"  pass 1: GT", flush=True)
        n_gt = split_to_buckets(GT, f"{td}/gt")
        print(f"    {n_gt:,} rows", flush=True)

        totals = collections.Counter()
        all_pval_log10_diffs = []
        all_score_diffs = []
        all_new_only_pv_log10 = []
        all_new_only_score = []
        for b in range(N_BUCKETS):
            r = compare_bucket(f"{td}/new.{b}.tsv", f"{td}/gt.{b}.tsv")
            for k, v in r.items():
                if isinstance(v, int):
                    totals[k] += v
            all_pval_log10_diffs.extend(r["pval_log10_diffs"])
            all_score_diffs.extend(r["score_diffs"])
            all_new_only_pv_log10.extend(r["new_only_pv_log10"])
            all_new_only_score.extend(r["new_only_score"])
            print(
                f"  bucket {b:2d}: new={r['new_hits']:>9,} gt={r['gt_hits']:>9,} "
                f"shared={r['shared']:>9,} new_only={r['new_only']:>8,}",
                flush=True,
            )

    sh = totals['shared']
    print()
    print(f"=== Totals over {N_BUCKETS} buckets (PER-HIT comparison; keyed by motif+alt+seq+start+stop+strand) ===")
    print(f"  NEW hits:                            {totals['new_hits']:>12,}")
    print(f"  GT hits:                             {totals['gt_hits']:>12,}")
    print(f"  shared (same position):              {totals['shared']:>12,}")
    print(f"  NEW-only:                            {totals['new_only']:>12,}")
    print(f"  GT-only:                             {totals['gt_only']:>12,}")
    print()
    print(f"=== On shared hits ({sh:,}) ===")
    if sh:
        print(f"  score byte-equal:                    {totals['score_eq_str']:>12,}  ({100*totals['score_eq_str']/sh:.4f}%)")
        print(f"  score equal within 1e-6 rel tol:     {totals['score_eq_close']:>12,}  ({100*totals['score_eq_close']/sh:.4f}%)")
        print(f"  p-value byte-equal:                  {totals['pval_eq_str']:>12,}  ({100*totals['pval_eq_str']/sh:.4f}%)")
        print(f"  p-value equal within 1e-6 rel tol:   {totals['pval_eq_close']:>12,}  ({100*totals['pval_eq_close']/sh:.4f}%)")

    # Histograms
    if all_pval_log10_diffs:
        bins = [-30, -20, -15, -10, -8, -7, -6, -5, -4, -3, -2, -1, 0, 5]
        histogram(all_pval_log10_diffs, bins, "log10 |Δp| distribution on shared-hit p-value mismatches")
    if all_score_diffs:
        bins = [0, 1e-9, 1e-6, 1e-4, 1e-3, 0.01, 0.1, 1, 10]
        # log-spaced; transform
        log_diffs = [math.log10(d) if d > 0 else -30 for d in all_score_diffs]
        bins_l = [-30, -15, -9, -6, -4, -3, -2, -1, 0, 5]
        histogram(log_diffs, bins_l, "log10 |Δscore| on shared-hit score mismatches")
    if all_new_only_pv_log10:
        bins = [-30, -20, -15, -10, -8, -7, -6, -5, -4.5, -4.2, -4.1, -4.0]
        histogram(all_new_only_pv_log10, bins, "log10 p-value of NEW-only hits (threshold was 1e-4 = -4)")
    if all_new_only_score:
        bins = [0, 5, 10, 15, 20, 25, 30, 50, 100, 1000]
        histogram(all_new_only_score, bins, "score distribution of NEW-only hits")


if __name__ == "__main__":
    main()
