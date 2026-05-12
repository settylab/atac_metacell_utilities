"""Compare gene_x_tf matrices produced by the downstream chain on two FIMO
inputs: FIMO 5.5.5 (this branch) vs FIMO 5.1.1 ground truth (2023).

Inputs (CSVs written by compute_gene_tf.py):
  results/gxtf-5p5p5/gene_x_tf/gene_x_tf.csv  -- from FIMO 5.5.5 output
  results/gxtf-5p1p1/gene_x_tf/gene_x_tf.csv  -- from 2023 GT FIMO 5.1.1

CSV layout: rows = genes (sc_rna.var_names), cols = TFs (sc_atac.uns['InSilicoChipColumns'])
Values: gene_x_peak (binarized gp_corrs) . peak_x_tf (varm['InSilicoChip']).

Reports:
  1. Shape diff:           genes/TFs gained/lost between the two matrices.
  2. Overall correlation:  Pearson + Spearman over common (gene, TF) entries.
  3. Per-TF column corr:   distribution of corr(col_5.5.5, col_5.1.1) across TFs,
                           weighted by nonzero entries.
  4. Top-N divergent TFs:  TFs whose column correlation < threshold or whose
                           column mean differs substantially.
  5. Scatter plots:        for a handful of top-divergent TFs, plot 5.5.5 vs 5.1.1
                           one dot per gene.

Outputs:
  exploration/fimo-version-investigation/figures/  (PNGs)
  Prints a markdown-format summary to stdout, also writes it to
    exploration/fimo-version-investigation/figures/summary.md
"""
from __future__ import annotations

import os
import sys
import math
import pathlib

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

REPO = pathlib.Path(__file__).resolve().parents[2]
A_PATH = REPO / "results" / "gxtf-5p5p5" / "gene_x_tf" / "gene_x_tf.csv"
B_PATH = REPO / "results" / "gxtf-5p1p1" / "gene_x_tf" / "gene_x_tf.csv"
OUT = REPO / "exploration" / "fimo-version-investigation" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

print("# gene_x_tf comparison: FIMO 5.5.5 (NEW) vs FIMO 5.1.1 GT")
print()
print(f"NEW: `{A_PATH.relative_to(REPO)}`")
print(f"GT:  `{B_PATH.relative_to(REPO)}`")
print()


def load_csv(path):
    print(f"reading {path} ...", file=sys.stderr)
    df = pd.read_csv(path, index_col=0)
    print(f"  shape: {df.shape} (genes x TFs)", file=sys.stderr)
    return df


df_a = load_csv(A_PATH)  # 5.5.5
df_b = load_csv(B_PATH)  # 5.1.1 GT

# 1) Shape diff
genes_a, genes_b = set(df_a.index), set(df_b.index)
tfs_a, tfs_b = set(df_a.columns), set(df_b.columns)
genes_common = sorted(genes_a & genes_b)
tfs_common = sorted(tfs_a & tfs_b)
print("## 1. Shape diff")
print()
print(f"|  | NEW (5.5.5) | GT (5.1.1) | common | NEW-only | GT-only |")
print(f"|---|---:|---:|---:|---:|---:|")
print(f"| genes | {len(genes_a):,} | {len(genes_b):,} | {len(genes_common):,} | {len(genes_a - genes_b):,} | {len(genes_b - genes_a):,} |")
print(f"| TFs   | {len(tfs_a):,} | {len(tfs_b):,} | {len(tfs_common):,} | {len(tfs_a - tfs_b):,} | {len(tfs_b - tfs_a):,} |")
print()

# Subset both to the common gene-TF grid for direct comparison
A = df_a.loc[genes_common, tfs_common].astype(np.float64)
B = df_b.loc[genes_common, tfs_common].astype(np.float64)
print(f"common subset: {A.shape[0]:,} genes x {A.shape[1]:,} TFs = "
      f"{A.size:,} entries each")
print()

# 2) Overall correlation
print("## 2. Overall correlation (flattened over common (gene, TF) entries)")
print()
a_flat = A.values.ravel()
b_flat = B.values.ravel()

nz_mask = (a_flat != 0) | (b_flat != 0)
nz_a = a_flat[nz_mask]
nz_b = b_flat[nz_mask]

pear_all, _ = stats.pearsonr(a_flat, b_flat)
pear_nz, _ = stats.pearsonr(nz_a, nz_b)
# spearman over full flattened is slow & dominated by zeros; restrict to nz
spear_nz, _ = stats.spearmanr(nz_a, nz_b)

n_entries = len(a_flat)
print(f"| Domain | n | Pearson | Spearman |")
print(f"|---|---:|---:|---:|")
print(f"| all entries  | {n_entries:,} | {pear_all:.4f} | (not computed, dominated by zeros) |")
print(f"| nonzero in either | {len(nz_a):,} | {pear_nz:.4f} | {spear_nz:.4f} |")
print()

# Zero-rate
zero_a = (a_flat == 0).mean()
zero_b = (b_flat == 0).mean()
print(f"zero-fraction NEW: {zero_a:.4%}  GT: {zero_b:.4%}")
print()

# 3) Per-TF column correlation distribution
print("## 3. Per-TF column correlation")
print()
print("For each TF, compute Pearson corr of the 2 columns (over common genes).")
print()

per_tf = []
for tf in tfs_common:
    a_col = A[tf].values
    b_col = B[tf].values
    nz = (a_col != 0) | (b_col != 0)
    if nz.sum() < 3:
        per_tf.append((tf, np.nan, int(nz.sum()), float(a_col.sum()), float(b_col.sum())))
        continue
    if a_col[nz].std() < 1e-12 or b_col[nz].std() < 1e-12:
        per_tf.append((tf, np.nan, int(nz.sum()), float(a_col.sum()), float(b_col.sum())))
        continue
    r, _ = stats.pearsonr(a_col[nz], b_col[nz])
    per_tf.append((tf, r, int(nz.sum()), float(a_col.sum()), float(b_col.sum())))

per_tf_df = pd.DataFrame(per_tf, columns=["tf", "corr", "nz", "sum_a", "sum_b"]).set_index("tf")
per_tf_df["mean_a"] = per_tf_df["sum_a"] / A.shape[0]
per_tf_df["mean_b"] = per_tf_df["sum_b"] / A.shape[0]
per_tf_df["mean_ratio"] = per_tf_df["sum_a"] / per_tf_df["sum_b"].replace(0, np.nan)
per_tf_df = per_tf_df.sort_values("corr")
per_tf_df.to_csv(OUT / "per_tf_correlations.csv")

valid = per_tf_df["corr"].dropna()
print(f"TFs scored: {len(valid):,} / {len(tfs_common):,}")
print(f"per-TF column corr: median={valid.median():.4f}  IQR=({valid.quantile(.25):.4f}, {valid.quantile(.75):.4f})  min={valid.min():.4f}  max={valid.max():.4f}")
print(f"TFs with corr < 0.80: {(valid < 0.80).sum():,}")
print(f"TFs with corr < 0.50: {(valid < 0.50).sum():,}")
print(f"TFs with corr < 0.20: {(valid < 0.20).sum():,}")
print()

# 4) Top-divergent TFs
print("## 4. Top-divergent TFs (lowest column correlation)")
print()
print("| TF | corr | nz | mean_5.5.5 | mean_5.1.1 | sum_ratio (5.5.5/5.1.1) |")
print("|---|---:|---:|---:|---:|---:|")
for tf, row in per_tf_df.head(20).iterrows():
    print(f"| `{tf}` | {row['corr']:.4f} | {int(row['nz']):,} | "
          f"{row['mean_a']:.4g} | {row['mean_b']:.4g} | "
          f"{row['mean_ratio']:.3g} |")
print()

# === Figures ===

# Figure 1: histogram of per-TF column correlations
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.hist(valid.values, bins=50, color="#4477AA", edgecolor="white")
ax.axvline(valid.median(), color="red", linestyle="--",
           label=f"median = {valid.median():.3f}")
ax.set_xlabel("per-TF Pearson corr (NEW 5.5.5 vs GT 5.1.1)")
ax.set_ylabel("number of TFs")
ax.set_title(f"Per-TF column correlation distribution (n={len(valid):,} TFs)")
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "per_tf_corr_hist.png", dpi=150)
plt.close(fig)
print(f"saved: {OUT / 'per_tf_corr_hist.png'}", file=sys.stderr)

# Figure 2: scatter of TF column means (5.5.5 vs 5.1.1)
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(per_tf_df["mean_b"], per_tf_df["mean_a"], s=6, alpha=0.4, color="#4477AA")
lim = max(per_tf_df[["mean_a", "mean_b"]].max().max(), 1e-9)
ax.plot([0, lim], [0, lim], "r--", lw=0.8, label="y = x")
ax.set_xlabel("TF column mean — GT (FIMO 5.1.1)")
ax.set_ylabel("TF column mean — NEW (FIMO 5.5.5)")
ax.set_title(f"TF column means (n={len(per_tf_df):,} common TFs)")
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "tf_column_mean_scatter.png", dpi=150)
plt.close(fig)
print(f"saved: {OUT / 'tf_column_mean_scatter.png'}", file=sys.stderr)

# Figure 3: scatter of top-3 most-divergent TFs (gene-level)
n_show = min(6, len(per_tf_df))
top_div = per_tf_df.dropna(subset=["corr"]).head(n_show).index.tolist()
fig, axes = plt.subplots(2, 3, figsize=(13, 8))
for ax, tf in zip(axes.ravel(), top_div):
    a_col = A[tf].values
    b_col = B[tf].values
    ax.scatter(b_col, a_col, s=5, alpha=0.4, color="#444")
    lim = max(abs(a_col).max(), abs(b_col).max(), 1e-9)
    ax.plot([-lim, lim], [-lim, lim], "r--", lw=0.5)
    r = per_tf_df.loc[tf, "corr"]
    nz = int(per_tf_df.loc[tf, "nz"])
    ax.set_title(f"{tf} | corr={r:.3f} | nz={nz}")
    ax.set_xlabel("GT (5.1.1)")
    ax.set_ylabel("NEW (5.5.5)")
for ax in axes.ravel()[len(top_div):]:
    ax.axis("off")
fig.suptitle("Top-divergent TFs: gene-level scatter (NEW vs GT)", y=1.02)
fig.tight_layout()
fig.savefig(OUT / "top_divergent_tf_scatter.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"saved: {OUT / 'top_divergent_tf_scatter.png'}", file=sys.stderr)

# Figure 4: overall density-style hexbin of NEW vs GT entries (nonzero in either)
fig, ax = plt.subplots(figsize=(6, 6))
hb = ax.hexbin(nz_b, nz_a, gridsize=50, bins="log", cmap="viridis", mincnt=1)
lim = max(abs(nz_a).max(), abs(nz_b).max(), 1e-9)
ax.plot([-lim, lim], [-lim, lim], "r--", lw=0.5, label="y = x")
ax.set_xlabel("GT (FIMO 5.1.1) entry")
ax.set_ylabel("NEW (FIMO 5.5.5) entry")
ax.set_title(f"All non-zero entries (n={len(nz_a):,})  Pearson={pear_nz:.4f}  Spearman={spear_nz:.4f}")
plt.colorbar(hb, ax=ax, label="log10(count)")
ax.legend(loc="upper left")
fig.tight_layout()
fig.savefig(OUT / "all_entries_hexbin.png", dpi=150)
plt.close(fig)
print(f"saved: {OUT / 'all_entries_hexbin.png'}", file=sys.stderr)

print()
print(f"all figures in: {OUT.relative_to(REPO)}")
