# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Snakemake pipeline for analyzing ATAC-seq and RNA-seq data at the metacell level using the SEACells algorithm. Produces:
- In-silico ChIP transcription factor targets
- ChromVAR accessibility results
- Primed and lineage-specific peaks in single-cell ATAC data

## Environment Setup

- **Python 3.8** via conda: `conda env create -n atac-metacell-utils --file envs/environment.yaml`
- **R 4.1.0** packages via renv: `snakemake --cores 1 renv_init_restore`
- **MEME Suite 5.1.1** required externally (provides `fimo` tool)
- Key Python deps: scanpy, anndata, SEACells (GitHub), Palantir (GitHub), snakemake 7.x
- Key R deps: edgeR, chromVAR (via renv.lock in `envs/`)

## Architecture

### Pipeline Flow (two main pathways)

**Peak-Scoring Pathway** (primed & lineage-specific peaks):
`prepare_peak_file` → `all_seqs` → `fimo` → `gp_corr` + `diff_acc` + `ct_open_peaks` → `peak_selection` → `peak_scores`

**In-Silico ChIP / ChromVAR Pathway:**
`prepare_peak_file` → `all_seqs` → `fimo` → `peak_tf` → `compute_ins_chip` → `prep_chromvar` → `chromvar`
`compute_ins_chip` + `gp_corr` → `gene_x_tf`

### Key Directories

- `workflow/Snakefile` — Main pipeline definition with all rule dependencies
- `workflow/scripts/` — Individual rule scripts (Python and R). Some rules use `script:` directive (access snakemake object directly), others use `shell:` with argparse
- `config/config.yaml` — Template config; project-specific configs (e.g., `config/D1_Config.yaml`) override it via `--configfile`
- `src/` — Utility modules (density analysis, plotting) used outside the main pipeline
- `envs/` — Environment specs (conda yaml, pip requirements, renv lock)
- `data/` — Reference files (MEME motifs, GTF genome annotations)

### Data Model

All pipeline outputs are stored as layers/annotations on AnnData (h5ad) objects following scanpy conventions:
- `.varm[]` for peak-level annotations (OpenPeaks, FIMO, InSilicoChip, diff_acc, geneXTF)
- `.layers[]` for cell×gene matrices (primed/lineage_specific scores)
- `.obsm[]` for cell-level embeddings (chromVAR_deviations)
- `.uns[]` for metadata (gp_corrs, column names)

### Configuration

Config requires four AnnData paths (`atac`, `rna`, `sc_atac`, `sc_rna`), output directory, genome, and cell type comparison parameters. The `diff_acc` and `peak_selection` config sections are dataset-specific and must be set per experiment. Config can be overridden via `--config key=value` or `--configfile`.

### Sentinel Files

Many rules use empty sentinel files (`.peak_tf`, `.gp_corrs`, `.diff_acc`, etc.) created with `touch` to signal completion, since actual outputs are directories with multiple files.
