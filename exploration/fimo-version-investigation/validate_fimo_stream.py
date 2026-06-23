"""Streaming variant of validate_fimo.py that works with limited memory.

Two-pass approach:
1. Pass 1 over GT: build a sorted-by-key list of (motif_id, motif_alt_id,
   sequence_name) -> (score, pvalue) and write to a temp file.
2. Pass 1 over NEW: same.
3. Sort both temp files by key and join via a streaming merge.

For 60M rows, this is feasible because each side fits 60M small tuples in
~6 GB RAM, but doing both sides at once + extra dicts blows past head-node
norms. So we hash-bucket each side by sequence_name into 32 buckets on
disk and compare bucket-by-bucket.
"""
import os, sys, hashlib, tempfile, collections

NEW = "/fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/results/tcell-bm-validate/fimo_result/fimo.tsv"
GT  = "/fh/fast/setty_m/user/yhuang2/cell_cell_communication/tcell_depleted_bone_marrow/data/fimo.tsv"
N_BUCKETS = 16

def bucket(name): return int(hashlib.md5(name.encode()).hexdigest()[:4], 16) % N_BUCKETS

def split_to_buckets(path, prefix):
    """Stream rows into N_BUCKETS files keyed by sequence_name hash.
    Each output line: motif_id\\tmotif_alt_id\\tsequence_name\\tscore\\tpvalue\\n"""
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
                # cols: motif_id, motif_alt_id, seq_name, start, stop, strand, score, pval, qval, matched
                if len(cols) < 8:
                    continue
                seq_name = cols[2]
                b = bucket(seq_name)
                handles[b].write(f"{cols[0]}\t{cols[1]}\t{seq_name}\t{cols[6]}\t{cols[7]}\n")
                n += 1
    finally:
        for h in handles: h.close()
    return n

def compare_bucket(new_path, gt_path):
    """For one bucket: build sets of (motif, alt, seq, score) and (motif, alt, seq, pval),
    return per-bucket counts."""
    def read(p):
        keys = set(); pvs = {}
        with open(p) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 5: continue
                keys.add((cols[0], cols[1], cols[2], cols[3]))
                pvs[(cols[0], cols[1], cols[2])] = cols[4]
        return keys, pvs
    new_k, new_p = read(new_path)
    gt_k, gt_p = read(gt_path)
    common = new_k & gt_k
    n_only = new_k - gt_k
    gt_only = gt_k - new_k
    pkeys = set(new_p) & set(gt_p)
    pmismatch = sum(1 for k in pkeys if new_p[k] != gt_p[k])
    return {
        "new_keys": len(new_k), "gt_keys": len(gt_k), "common_score_keys": len(common),
        "new_only": len(n_only), "gt_only": len(gt_only),
        "pval_shared": len(pkeys), "pval_mismatch": pmismatch,
    }

scratch = os.environ.get("SCRATCH_DIR", "/tmp")
with tempfile.TemporaryDirectory(prefix="fimo_diff_", dir=scratch) as td:
    print(f"bucketing into {td}/", flush=True)
    print(f"  pass 1: NEW", flush=True)
    n_new = split_to_buckets(NEW, f"{td}/new")
    print(f"    {n_new:,} rows", flush=True)
    print(f"  pass 1: GT", flush=True)
    n_gt = split_to_buckets(GT, f"{td}/gt")
    print(f"    {n_gt:,} rows", flush=True)

    totals = collections.Counter()
    for b in range(N_BUCKETS):
        r = compare_bucket(f"{td}/new.{b}.tsv", f"{td}/gt.{b}.tsv")
        for k, v in r.items(): totals[k] += v
        print(f"  bucket {b:2d}: new={r['new_keys']:,} gt={r['gt_keys']:,} common={r['common_score_keys']:,}", flush=True)

    print()
    print(f"=== Totals over {N_BUCKETS} buckets ===")
    for k, v in totals.items():
        print(f"  {k}: {v:,}")
    if totals['gt_keys']:
        ms = totals['common_score_keys'] / totals['gt_keys']
        print(f"\n  (motif_id, motif_alt_id, sequence_name, score) match rate: {totals['common_score_keys']:,} / {totals['gt_keys']:,} ({100*ms:.2f}%)")
    if totals['pval_shared']:
        pmatch = totals['pval_shared'] - totals['pval_mismatch']
        print(f"  p-value match rate (on (motif,alt,seq) shared keys): {pmatch:,} / {totals['pval_shared']:,} ({100*pmatch/totals['pval_shared']:.4f}%)")
