"""
Differential accessibility of metacell (SEACell) pseudo-bulk ATAC data.

This single script replaces the previous three-stage pipeline:

    meta_mtx.py  ->  diff_acc.R  ->  diff_acc.py

It keeps the snakemake *inputs* of meta_mtx.py (the single-cell and metacell
AnnData objects pulled from ``config['anndata']``, the SEACell label, the
cell-type grouping column, and the modality), and the snakemake *outputs* of
diff_acc.py (the ``<groupA>_<groupB>_diff_acc`` tables written into
``atac.varm`` of the metacell AnnData, plus the sentinel ``output[0]``). It also
reproduces the per-comparison MA-plot PDFs that diff_acc.R wrote alongside its
results (``<groupA>_<groupB>_MA.pdf``).

The intermediate ``.mtx``/``.csv`` round-trip that used to ferry data from
Python into R is gone: the edgeR GLM quasi-likelihood workflow is run directly
in-process with edgePython (https://github.com/pachterlab/edgePython), so the
counts never leave memory.

Snakemake interface
-------------------
config['anndata']['sc_atac']     single-cell ATAC AnnData path
config['anndata']['atac']        metacell ATAC AnnData path (read + written back)
config['diff_acc']['cell_type_obs']  .obs column used for grouping (group_variable)
config['modality']               modality label (informational)
params.seacell_label             .obs column holding the SEACell assignment
params.to_compare                "A,B/C,D" -> compare (A vs B) and (C vs D)
params.out_dir                   (optional) directory for the MA-plot PDFs;
                                 defaults to the directory of output[0]
output[0]                        sentinel file touched on completion
"""

import os

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import edgepython as ep

import matplotlib

matplotlib.use("Agg")  # headless backend for cluster / batch use
import matplotlib.pyplot as plt


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def signif(values, digits=3):
    """Round to ``digits`` significant figures (R's ``signif``)."""
    x = np.asarray(values, dtype=float)
    out = np.array(x, dtype=float)
    nz = np.isfinite(x) & (x != 0)
    if np.any(nz):
        mags = np.floor(np.log10(np.abs(x[nz])))
        factor = 10.0 ** (digits - 1 - mags)
        out[nz] = np.round(x[nz] * factor) / factor
    return out


def plot_ma(logcpm, logfc, out, pair, path, lfc_cutoff=1.25, padj_cutoff=0.01):
    """Write an MA plot (logFC vs logCPM) for one comparison.

    Faithful port of the plotting block in diff_acc.R: all tested peaks in
    grey, group-A peaks (logFC <= -lfc_cutoff & padj <= padj_cutoff) in red,
    group-B peaks (logFC >= lfc_cutoff & padj <= padj_cutoff) in blue. Point
    coordinates use the unrounded logCPM/logFC of the tested peaks, while
    significance membership is taken from the rounded output table, matching R.
    """
    group_a, group_b = pair[0], pair[1]

    sig_a = out.index[(out["logFC"] <= -lfc_cutoff) & (out["padj_fdr"] <= padj_cutoff)]
    sig_b = out.index[(out["logFC"] >= lfc_cutoff) & (out["padj_fdr"] <= padj_cutoff)]
    sig_a = [f for f in sig_a if f in logcpm.index]
    sig_b = [f for f in sig_b if f in logcpm.index]

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(logcpm.values, logfc.values, s=2, c="lightgray", linewidths=0)
    if sig_a:
        ax.scatter(
            logcpm[sig_a].values, logfc[sig_a].values,
            s=12, c="red", linewidths=0, label=f"{group_a} peaks",
        )
    if sig_b:
        ax.scatter(
            logcpm[sig_b].values, logfc[sig_b].values,
            s=12, c="dodgerblue", linewidths=0, label=f"{group_b} peaks",
        )
    ax.axhline(0, color="black", linewidth=0.5, zorder=0)
    ax.set_xlabel("logCPM")
    ax.set_ylabel("logFC")
    ax.set_title(f"{group_a} vs {group_b}")
    if sig_a or sig_b:
        ax.legend(loc="lower left", frameon=False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def assign_metacell_celltypes(sc_ad, meta_ad, seacell_label, group_variable):
    """Label each metacell with the dominant cell type of its member cells.

    Equivalent to the original meta_mtx.py idiom
    ``obs[group].groupby(seacell).value_counts().groupby(level=0).head(1)``
    but written to be robust across pandas versions.
    """
    member = pd.DataFrame(
        {
            "seacell": np.asarray(sc_ad.obs[seacell_label].values),
            "celltype": np.asarray(sc_ad.obs[group_variable].values),
        }
    )
    # Dominant (most frequent) cell type per SEACell; ties resolved by the
    # order value_counts returns, matching the original .head(1) behaviour.
    dominant = member.groupby("seacell")["celltype"].agg(
        lambda s: s.value_counts().index[0]
    )
    meta_ad.obs[group_variable] = (
        dominant.reindex(meta_ad.obs_names).astype(str).values
    )
    return meta_ad


def diff_acc_edgepython(
    sub,
    peak_names,
    groups,
    pair,
    detection_thresh=0.3,
    lfc_round=3,
    fdr_signif=3,
    mean_round=2,
    ma_plot_path=None,
    lfc_cutoff=1.25,
    padj_cutoff=0.01,
):
    """Run the edgeR QL differential-accessibility test for one comparison.

    Faithful port of ``diff_acc_edgeR`` from diff_acc.R.

    Parameters
    ----------
    sub : ndarray (metacells x peaks)
        Dense raw counts for the metacells belonging to either group of ``pair``.
    peak_names : array-like
        Peak (feature) names, aligned to the columns of ``sub``.
    groups : ndarray
        Group label for each row (metacell) of ``sub``.
    pair : sequence of str
        ``[groupA, groupB]``. groupA is the reference level, so a positive
        ``logFC`` means higher accessibility in groupB.
    ma_plot_path : str, optional
        If given, an MA-plot PDF is written to this path.
    lfc_cutoff, padj_cutoff : float
        |logFC| and FDR thresholds used to colour significant peaks in the MA
        plot (R defaults: 1.25 and 0.01).

    Returns
    -------
    pandas.DataFrame indexed by feature, with columns
    ``feature, logFC, padj_fdr, groupA_N, groupB_N, mean_groupA, mean_groupB``
    covering *all* peaks (untested peaks get logFC=0, padj_fdr=1).
    """
    group_a, group_b = pair[0], pair[1]
    sub = np.asarray(sub, dtype=np.float64)
    peak_names = np.asarray(peak_names)

    mask_a = groups == group_a
    mask_b = groups == group_b
    n_a, n_b = int(mask_a.sum()), int(mask_b.sum())
    if n_a < 1 or n_b < 1:
        raise ValueError(
            f"Comparison {group_a} vs {group_b}: missing metacells "
            f"(n_{group_a}={n_a}, n_{group_b}={n_b})."
        )

    # Per-metacell library size over ALL peaks (R: colSums of the subset).
    libsize = sub.sum(axis=1, keepdims=True)
    libsize[libsize == 0] = 1.0  # guard against empty metacells
    # logCPM-style values, used only for the reported per-group means.
    logcounts = np.log(1e6 * sub / libsize + 1.0)

    mean_a = logcounts[mask_a].mean(axis=0)
    mean_b = logcounts[mask_b].mean(axis=0)

    # Detection rate per group, rounded to 2 dp before thresholding (as in R).
    det_a = np.round((sub[mask_a] > 0).mean(axis=0), 2)
    det_b = np.round((sub[mask_b] > 0).mean(axis=0), 2)
    feature_mask = (det_a >= detection_thresh) | (det_b >= detection_thresh)

    if not feature_mask.any():
        raise ValueError(
            f"Comparison {group_a} vs {group_b}: no peaks pass the "
            f"detection-rate filter (>= {detection_thresh})."
        )

    # Cellular detection rate per metacell over ALL peaks -> design covariate.
    cdr = (sub > 0).mean(axis=1)

    # --- edgeR / edgePython GLM quasi-likelihood pipeline --------------------
    counts_filt = sub[:, feature_mask].T          # peaks_filt x metacells
    feat_names_filt = peak_names[feature_mask]

    # Design ~ cdr + group, reference level = group_a (so the group
    # coefficient -- the last column -- is the B-vs-A effect that is tested).
    group_dummy = mask_b.astype(float)
    design = np.column_stack([np.ones(sub.shape[0]), cdr, group_dummy])

    # No calc_norm_factors() on purpose: the original R kept norm.factors = 1
    # and lib.size = colSums (the plain DGEList default), which edgePython
    # reproduces. Add ep.calc_norm_factors() here if TMM is ever desired.
    y = ep.make_dgelist(counts=counts_filt)
    y = ep.estimate_disp(y, design)
    fit = ep.glm_ql_fit(y, design)
    res = ep.glm_ql_ftest(fit)                    # default coef = last = group
    table = ep.top_tags(res, n=counts_filt.shape[0])["table"]

    # top_tags indexes rows by their position in counts_filt; map back to names.
    tested_names = feat_names_filt[np.asarray(table.index, dtype=int)]
    logfc = pd.Series(
        np.round(table["logFC"].to_numpy(), lfc_round), index=tested_names
    )
    padj = pd.Series(signif(table["FDR"].to_numpy(), fdr_signif), index=tested_names)

    # Unrounded coordinates kept for the MA plot (R plots lrt$table directly).
    logcpm_tested = pd.Series(table["logCPM"].to_numpy(), index=tested_names)
    logfc_tested = pd.Series(table["logFC"].to_numpy(), index=tested_names)

    # --- Assemble one row per peak (R merges with all.y = TRUE) --------------
    peak_index = pd.Index(peak_names)
    out = pd.DataFrame(index=peak_index)
    out["feature"] = peak_names
    out["logFC"] = logfc.reindex(peak_index).fillna(0.0).to_numpy()
    out["padj_fdr"] = padj.reindex(peak_index).fillna(1.0).to_numpy()
    out["groupA_N"] = n_a
    out["groupB_N"] = n_b
    out["mean_groupA"] = np.round(mean_a, mean_round)
    out["mean_groupB"] = np.round(mean_b, mean_round)

    out = out.sort_values("padj_fdr", kind="mergesort")
    out.index = out["feature"]
    out.index.name = None

    if ma_plot_path is not None:
        plot_ma(
            logcpm_tested, logfc_tested, out, pair, ma_plot_path,
            lfc_cutoff=lfc_cutoff, padj_cutoff=padj_cutoff,
        )

    return out


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #
def main():
    cfg = snakemake.config

    sc_atac_path = cfg["anndata"]["sc_atac"]
    atac_path = cfg["anndata"]["atac"]
    group_variable = cfg["diff_acc"]["cell_type_obs"]
    modality = cfg.get("modality", "atac") if hasattr(cfg, "get") else cfg["modality"]
    seacell_label = snakemake.params["seacell_label"]

    # "A,B/C,D" -> [["A", "B"], ["C", "D"]]
    comparisons = [
        pair.split(",") for pair in str(snakemake.params["to_compare"]).split("/")
    ]

    # Sentinel output path (formerly diff_acc.py output[0]).
    sentinel_path = str(snakemake.output[0])

    # MA-plot directory: prefer an explicit out_dir param (mirroring the
    # original diff_acc.py, which read params["out_dir"][0]); otherwise fall
    # back to the directory holding the sentinel output.
    try:
        out_dir = snakemake.params["out_dir"]
        ma_dir = out_dir[0] if isinstance(out_dir, (list, tuple)) else out_dir
    except (KeyError, TypeError, IndexError):
        ma_dir = os.path.dirname(sentinel_path)
    ma_dir = ma_dir or "."
    os.makedirs(ma_dir, exist_ok=True)

    print("Loading data...")
    sc_ad = sc.read(sc_atac_path)
    meta_ad = sc.read(atac_path)

    # Dominant cell type per metacell (formerly meta_mtx.py).
    meta_ad = assign_metacell_celltypes(
        sc_ad, meta_ad, seacell_label, group_variable
    )

    raw = meta_ad.layers["raw"]
    peak_names = meta_ad.var_names.to_numpy()
    groups_all = meta_ad.obs[group_variable].to_numpy()

    print(f"Performing differential accessibility ({modality}) with edgePython...")
    for pair in comparisons:
        group_a, group_b = pair[0], pair[1]
        print(f"  {group_a} vs {group_b}")

        sel = np.isin(groups_all, [group_a, group_b])
        sub = raw[sel]
        sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)

        ma_plot_path = os.path.join(ma_dir, f"{group_a}_{group_b}_MA.pdf")
        result = diff_acc_edgepython(
            sub, peak_names, groups_all[sel], pair, ma_plot_path=ma_plot_path
        )

        # Store aligned to var_names, exactly as the old diff_acc.py did.
        meta_ad.varm[f"{group_a}_{group_b}_diff_acc"] = result.loc[
            meta_ad.var_names, :
        ]

    print("Writing updated AnnData...")
    meta_ad.write(atac_path)

    # Touch the sentinel output.
    os.makedirs(os.path.dirname(sentinel_path) or ".", exist_ok=True)
    with open(sentinel_path, "a"):
        os.utime(sentinel_path, None)


if __name__ == "__main__":
    main()
