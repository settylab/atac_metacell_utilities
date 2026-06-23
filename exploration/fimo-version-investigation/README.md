# FIMO version investigation (May 2026)

Scripts and figures from the validation-and-comparison work on issue
[settylab/sarah-nexus#4](https://github.com/settylab/sarah-nexus/issues/4):
porting `atac_metacell_utilities` to a pure-Python pipeline, validating
the new FIMO 5.5.5 output against the 2023 FIMO 5.1.1 ground truth (GT),
and then quantifying the downstream effect on `gene_x_tf`.

These scripts were originally written under `results/tcell-bm-validate/`
and moved here for tidiness. Data references inside the scripts still
point to that directory (the data itself was not moved); the wrappers
that call each Python script were updated to point at the new
`exploration/fimo-version-investigation/` location of the `.py` files.

## Layout

```
exploration/fimo-version-investigation/
├── README.md                       <- this file
├── compare_gene_x_tf.py            <- TASK 1: gene_x_tf NEW vs GT comparison
├── figures/                        <- all PNGs from compare_gene_x_tf.py
│   ├── per_tf_corr_hist.png
│   ├── tf_column_mean_scatter.png
│   ├── top_divergent_tf_scatter.png
│   ├── all_entries_hexbin.png
│   ├── per_tf_correlations.csv     <- per-TF corr table sorted ascending
│   └── summary.md                  <- markdown summary printed by compare_gene_x_tf.py
└── (FIMO-version validation, prior session) ----------------------------
    ├── validate_fimo_perhit.py     <- per-hit byte-equal diff NEW vs GT,
    │                                  keyed by full (motif, alt, seq, start, stop, strand)
    ├── validate_fimo_stream.py     <- earlier "by triple" variant (methodology bug)
    ├── diff_chunk0_three_way.py    <- A (5.5.5) vs B (5.1.1) vs GT on chunk_0
    ├── threshold_what_if.py        <- what-if: tightening --thresh would cut GT coverage
    ├── seq_gl_for_diff.R           <- R seq_gl variant used during the byte-identical check
    ├── run_*.sh                    <- srun wrappers for the .py and .R scripts above
    └── (see "Script-by-script notes" below for the wrapper roles)
```

## Workflow

### TASK 1 — gene_x_tf comparison (compare_gene_x_tf.py)

Quantifies the downstream effect of swapping the FIMO 5.5.5 output
(from the new scatter-gather pipeline) for the 2023 FIMO 5.1.1 GT
output.

The two downstream chains were run in `results/gxtf-5p5p5/` and
`results/gxtf-5p1p1/`. Each contains its own copies of `sc_atac.h5ad`
and `sc_rna.h5ad` (since `peak_tf`, `compute_ins_chip`, and
`compute_gene_tf` all write back into those AnnDatas in place); the
meta-cell AnnDatas (`atac_mc`, `rna_mc`) are read-only and symlinked
from the upstream location.

Run:

```bash
srun --jobid=<your-allocation> --overlap --ntasks=1 --cpus-per-task=2 --mem=16G \
    bash -c '
    eval "$($MAMBA_EXE shell hook --shell bash --root-prefix $MAMBA_ROOT_PREFIX 2>/dev/null)"
    micromamba activate seacells
    python exploration/fimo-version-investigation/compare_gene_x_tf.py > exploration/fimo-version-investigation/figures/summary.md
    '
```

### FIMO-version validation (prior session)

These ran during the original port (May 8–11) and established that
the new FIMO 5.5.5 output **strictly contains** the 2023 GT byte-equal:
100% of GT hits are reproduced with the same score and p-value; 0 are
missing; the only difference is 9% of NEW-only hits at the `--thresh
1e-4` boundary. See issue #4 and the corresponding workspace reports
for the full audit trail.

## Script-by-script notes

`validate_fimo_perhit.py` — Two-pass, hash-bucketed streaming diff
keyed by `(motif_id, motif_alt_id, sequence_name, start, stop,
strand)`. Outputs counts and per-bucket / aggregated match rates for
score and p-value. Run via `run_perhit_diff.sh`. Replaces the earlier
`validate_fimo_stream.py` (which was keyed by `(motif, alt, seq)`
3-tuple — a methodology bug since multi-hit peaks would overwrite
p-values).

`diff_chunk0_three_way.py` — Restricts to chunk_0's peaks and compares
three datasets: A = new 5.5.5 chunk_0, B = 5.1.1 controlled re-run on
the same chunk_0 input, GT = 2023 GT filtered to chunk_0 peaks.
Established that the 2023 invocation used unrecorded extra flags
(probably `--max-stored-scores`) that today's 5.1.1 can't reproduce.
Run via `run_chunk0_diff.sh`.

`threshold_what_if.py` — Counter-factual analysis: would tightening
`--thresh 1e-4` to `1e-5` or `5e-5` recover the 26%-vs-100% gap with
the 2023 GT? Answer: no — tightening would catastrophically destroy
GT coverage (80% loss at 1e-5). Run via `run_threshold_whatif.sh`.

`run_bg_test.sh` — Probes whether per-chunk background drift (FIMO
auto-computes a background from each chunk's FASTA) explains the
NEW-only boundary hits. Spoiler: it does not — switching to a global
`all_seqs.bg` background changes hit counts marginally.

`run_full_fimo.sh` — Convenience wrapper that runs the full FIMO
scatter via Snakemake against the in-repo `config/config.yaml`.

`run_seq_gl_r.sh` — Runs the R seq_gl variant on the same peaks input,
producing `all_seqs.r.fa` for byte-identical comparison with the
Python port output (`all_seqs.fa`).

`run_fimo_probe.sh`, `run_fimo511_probe.sh`, `run_fimo511_noflags.sh`
— Per-chunk FIMO probes at 5.5.5 (default scatter flags), 5.1.1
(`--no-qvalue`), and 5.1.1 (no flags, matching the 2023 trailer
comment) respectively.

`run_validate_diff.sh` — Legacy wrapper for the earlier
`validate_fimo_stream.py` 3-tuple diff. Kept for record.

`run_perhit_diff.sh`, `run_chunk0_diff.sh`, `run_threshold_whatif.sh`
— srun wrappers for the Python scripts named above.
