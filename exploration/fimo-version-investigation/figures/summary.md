reading /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/results/gxtf-5p5p5/gene_x_tf/gene_x_tf.csv ...
  shape: (17226, 701) (genes x TFs)
reading /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/results/gxtf-5p1p1/gene_x_tf/gene_x_tf.csv ...
  shape: (17226, 701) (genes x TFs)
saved: /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/exploration/fimo-version-investigation/figures/per_tf_corr_hist.png
saved: /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/exploration/fimo-version-investigation/figures/tf_column_mean_scatter.png
saved: /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/exploration/fimo-version-investigation/figures/top_divergent_tf_scatter.png
saved: /fh/fast/setty_m/user/yhuang2/sarah-nexus/repos/atac_metacell_utilities_s/exploration/fimo-version-investigation/figures/all_entries_hexbin.png
# gene_x_tf comparison: FIMO 5.5.5 (NEW) vs FIMO 5.1.1 GT

NEW: `results/gxtf-5p5p5/gene_x_tf/gene_x_tf.csv`
GT:  `results/gxtf-5p1p1/gene_x_tf/gene_x_tf.csv`

## 1. Shape diff

|  | NEW (5.5.5) | GT (5.1.1) | common | NEW-only | GT-only |
|---|---:|---:|---:|---:|---:|
| genes | 17,226 | 17,226 | 17,226 | 0 | 0 |
| TFs   | 701 | 701 | 701 | 0 | 0 |

common subset: 17,226 genes x 701 TFs = 12,075,426 entries each

## 2. Overall correlation (flattened over common (gene, TF) entries)

| Domain | n | Pearson | Spearman |
|---|---:|---:|---:|
| all entries  | 12,075,426 | 0.9967 | (not computed, dominated by zeros) |
| nonzero in either | 2,376,170 | 0.9958 | 0.9887 |

zero-fraction NEW: 80.3223%  GT: 80.6924%

## 3. Per-TF column correlation

For each TF, compute Pearson corr of the 2 columns (over common genes).

TFs scored: 701 / 701
per-TF column corr: median=1.0000  IQR=(1.0000, 1.0000)  min=0.7062  max=1.0000
TFs with corr < 0.80: 2
TFs with corr < 0.50: 0
TFs with corr < 0.20: 0

## 4. Top-divergent TFs (lowest column correlation)

| TF | corr | nz | mean_5.5.5 | mean_5.1.1 | sum_ratio (5.5.5/5.1.1) |
|---|---:|---:|---:|---:|---:|
| `ZNF613` | 0.7062 | 3,086 | 0.003244 | 0.004461 | 0.727 |
| `ZNF182` | 0.7896 | 2,388 | 0.001842 | 0.0008527 | 2.16 |
| `ZNF432` | 0.9126 | 7,101 | 0.004077 | 0.00191 | 2.13 |
| `VEZF1` | 0.9186 | 6,671 | 0.003773 | 0.00237 | 1.59 |
| `ZNF529` | 0.9257 | 6,643 | 0.01149 | 0.00658 | 1.75 |
| `ZNF180` | 0.9362 | 6,660 | 0.009754 | 0.004649 | 2.1 |
| `ZNF548` | 0.9390 | 6,659 | 0.007884 | 0.004435 | 1.78 |
| `ZNF283` | 0.9477 | 6,919 | 0.007316 | 0.004431 | 1.65 |
| `ZNF267` | 0.9518 | 6,343 | 0.009299 | 0.006216 | 1.5 |
| `ZNF383` | 0.9538 | 5,946 | 0.007328 | 0.005211 | 1.41 |
| `IRF3` | 0.9550 | 6,699 | 0.002719 | 0.002128 | 1.28 |
| `THAP1` | 0.9557 | 7,025 | 0.003556 | 0.003139 | 1.13 |
| `ZNF519` | 0.9560 | 7,274 | 0.008243 | 0.005286 | 1.56 |
| `ZNF467` | 0.9568 | 7,564 | 0.0107 | 0.00599 | 1.79 |
| `ZNF571` | 0.9576 | 6,109 | 0.005772 | 0.004025 | 1.43 |
| `ZNF770` | 0.9577 | 9,052 | 0.007546 | 0.005447 | 1.39 |
| `ZBTB17` | 0.9601 | 8,100 | 0.005298 | 0.004309 | 1.23 |
| `ZNF880` | 0.9642 | 7,320 | 0.002413 | 0.001882 | 1.28 |
| `ZNF304` | 0.9697 | 7,605 | 0.006269 | 0.004427 | 1.42 |
| `ZNF235` | 0.9713 | 2,910 | 0.002353 | 0.001608 | 1.46 |


all figures in: exploration/fimo-version-investigation/figures
