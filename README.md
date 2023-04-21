# geneTF

This repository contains a [`snakemake`](https://snakemake.readthedocs.io/en/stable/) pipeline to produce `gene x TF` matrices representing the regulatory potential of TFs on genes, improved chromVAR deviations scores using custom MOTIF annotations, peak accessibility measurements, and more. 

A `snakemake` pipeline is constructed using a set of *rules* with defined inputs and outputs.

## Overview:

Below is a DAG showing the rule dependencies in the `snakemake` pipeline. View the [`Snakefile`](https://github.com/settylab/geneTF/blob/main/workflow/Snakefile) to see details on each rule.

![DAG of workflow](./dag.png)

## Types of `gene x TF` matrices produced

Three versions of the `gene x TF` matrix are produced, each with a different proxy for quantifying regulation:

* sum of gene accessibility based on ATAC peak counts
* sum of FIMO scores
* sum of gene accessibility, peaks weighted by FIMO scores

## Environment Setup
In order to run this pipeline, certain python and R packages need to be installed. Importantly, the `environment.yaml`/`requirements.txt` files include the `snakemake` package.

### Modules
We also rely on some packages through `module load` system provide by SciComp. Specifically:

1. `Anaconda3/2020.02`
    * This is just to allow the use of `conda`, free to use whatever version you like
2. `R/4.1.0-foss-2020b`
    * This is the version of R used throughout the pipeline. If you wish to use a different R version / installation, modify the `renv.lock` file and the following rules:
        * `renv_install`
        * `renv_init_restore`
        * `all_seqs`
    * The `renv.lock` file will also prove to have problems if you change the R-version...
3. `MEME/5.1.1-foss-2019b-Perl-5.30.0-Python-3.7.4`
    * This includes the `fimo` tool used to find MOTIFs. If you wish to use your own installation of the MEME Suite, modifcy the following rule:
        * `fimo`

### Python Modules
1. Load Anaconda:

```
module load Anaconda3/2020.02
```

2. Create the conda environment 

There are two options, using `pip` vs using `conda`

With `conda`:

```
envName=gene-TF

conda env create -n "$envName" --file envs/environment.yaml

conda activate "$envName"
```

With `pip`:

In an activated, bare `conda` environment:

```
pip install -r envs/requirements.txt
```

### R modules

To ensure the project runs smoothly, an `renv.lock` file is included under the `envs/` directory. Using the [`renv` package](https://rstudio.github.io/renv/index.html) makes it easy to distribute the required packages and their correct versions.

If `renv` is not yet installed, run:

```
snakemake --cores 1 renv_install --config renv_loc=<path/to/user/R/library>
```
**NOTE:** `--cores` can tell `snakemake` how many cores are available for use. It will automatically run the jobs in parallel, if possible.


This will run the `snakemake` rule that will:
1. Load the R module
2. Install the `renv` package into the library location specified above

Next, we can install all the required R packages:

```
snakemake --cores 1 renv_init_restore
```

## Modify the `config.yaml` file

For each project using this pipeline, parameters can be set in the configuration file located at:
```
config/config.yaml
```

In order for the pipeline to run, you will have to set the following parameters:

* `atac` *str* : Path to the SEACell-level ATAC AnnData with normalized and log-transformed counts
* `rna` *str* : Path to the SEACell-level RNA AnnData with normalized and log-transformed counts
* `sc_atac` *str* : Path to the single-cell-level ATAC AnnData with matching peaks to `atac` AnnData. Must have `nFrags` annotation in `.obs`
* `sc_ad` *str* : Path to the single-cell-level ATAC/RNA AnnData with low-dim embedding and relevant `.obs` annotations for `plot_tfs` rule

**NOTE:** The RNA and ATAC SEACell AnnDatas should have matched/common `obs_names`

*For plotting TFs after running chromVAR:*
* `tf_list` *list/str* : space-separated list of TFs to plot, no quotations
* `obs_annos` *list/str* : space-separated list of `sc_ad.obs` annotations to plot, no quotations

### Optional Parameters:
For the following parameters, default values can be used

*Peak and genome params:*
**NOTE:** Make sure the `in_meme` file is appropriate for the chosen `genome`
* `in_meme` *str* : Path to the `.meme` file for FIMO, `default=data/cis-bp-tf-information.meme`
    * The default meme file was created using the script: `workflow/scripts/make_meme.py`
    * Other motif databases, including `CIS-BP_2.00/Homo_sapiens.meme`, can be found at: `/fh/fast/setty_m/grp/motif_databases/`
* `width` *int* : Number of base pairs to resize ATAC peaks to. The summit will be the midpoint, `default=150`
* `genome` *str* : Genome to use, will be referenced in the `all_seqs`, `chromvar`, and `gp_corr` rules, `default=hg38`
    * Currently supported genomes include:  
        * hg38
        * hg19
        * hg18
        * mm9
        * mm10

*In silico ChIP params:*
* `verbose` *str/flag* : Whether to print info for each TF, `default=""` (False), set to `"--verbose"` for True
* `min_chip_score` *float* : minimum in silico ChIP score to pass filtering, `default=0.10`
* `min_peak_hits` *int* : minimum TFs passing in silico ChIP filters, `default=30`
* `vmax` *int* : vmax for plotting, vmin = -vmax, `default=3`
* `embedding` *str* : `sc_ad.obsm` field for plotting, `default=umap`

*Gene-Peak correlation params:*
* `n_jobs` *int* : Number of cores for the `gp_rule`, `default=1`
* `test_set` *str/flag* : Whether to generate a test set, `default=""` (False), set to `"--test_set"` for True 
* `n_genes` *int* : number of genes to use in test set, if applicable, `default=20`
* `min_corr` *float* : minimum gene-peak correlation to pass filtering, `default=0.0`
* `max_pval` *float* : max p-value for gene-peak correlation to pass filtering, `default=0.1`
* `min_peaks` *int* : minimum number of significant peaks to pass gene filtering,  `default=2`

*Accessibility metrics params:*
* `ld_embedding` *str* : `atac.obsm` field for nearest neighbor computation, `default=X_svd`
* `read_length` *int* : Fragment length, `default=147`
* `op_nbors` *int* : Number of neighbors for kNN, `default=3`
* `open_peak_pval` *float* : Nominal p-value cutoff for open peaks, `default=0.01`


### Alternative config specification

Instead of altering the `config.yaml` file directly, you can also pass the params using the keys in the config file on the command line.

For example if we wanted to generate a test_set:

```
snakemake --cores 1 all --config test_set="--test_set"
```

#### Passing in a different config file
In the Snakefile, the global variable, `config`, is set with the line:

```
configfile: "config/config.yaml"
```

The keys defined in the default `config/config.yaml` can be **overwritten** by passing a different config file (still `.yaml`). For example:

```
snakemake --cores 1 all --configfile <new_config_file.yaml> 
```

**NOTE:** All keys defined by the `configfile:` statement, the `--configfile` and/or `--config` command line arguments are part of the final `config` dictionary. 

*If two methods define the same key, command line arguments **overwrite** keys from the `configfile:` statement*

Snakemake will output a statement to make it clear to the user:

> Config file config/config.yaml is extended by additional config specified via the command line


## Optional: Download `hg38.gtf`

For convenience, I have included a rule to download the *hg38.gtf* file into a data directory for use in the `gp_corr` rule.

```
snakemake --cores 1 dl_hg38
```

If a different genome is being used, create a `data/` directory and place the `.gtf` file inside. Modify the `config.yaml`, respectively.

## Running the pipeline

After the environment has been set up and the configuration file is set, the pipeline is ready to run!

**NOTE:** Reccommended to run on `gizmo`, not `rhino`! Cluster integration with `snakemake` directly for this pipeline is in the works!

We can always conduct a "dry run" to test that all required inputs/parameters have been set before running the pipeline for real:

```
snakemake -nr all
```
Here:
 * `-n` flag tells `snakemake` to make a "dry run"
 * `-r` flag will ouput the reason each rule is running

To generate all files, run:

```
snakemake --cores 1 all
```

**NOTE:** In the case that a rule fails, it is not necessary to rerun rules that successfully completed. 

In fact, calling the `all` rule again will only run the rules for which their output is missing.
```
snakemake --cores 1 all
```

For more control, individual rules can also be called:

If you do this, make sure all the input dependencies have already been generated! `-n` can be appended for a dry run to determine if inputs are available.

```
snakemake --cores 1 name_of_rule
```

After modifications of inputs or parameters, you can force a rule to rerun using the `-R` flag:
```
snakemake --cores 1 -R name_of_rule
```

For even more control, the `--allowed-rules` argument can be used to pass a space-separated list of rules that are allowed to run. This is helpful when trying to avoid re-running time-intensive rules such as `fimo` and `gp_corr` when calling `all` or individual rules that get inputs from upstream rules. 

For example, the following will run/re-run the `chromvar` rule any of the allowed rules if necessary. This is useful if you have the `fimo` matrices already computed and do not want to run that rule.
```
snakemake --cores 1 chromvar --allowed_rules peak_tf compute_ins_chip prep_chromvar
```


### Submitting to Slurm
We can submit jobs to the cluster non-interactively by creating a shell script calling `snakemake`. An template script is shown below

Output (stdout and stderr) captured by Slurm is written to a `.log` file, specified with the `--output` option.

Here:

* `$1` : name of `snakemake` rule to run
* `$2` : number of cores to use

```
#!/bin/bash

#SBATCH --cpus-per-task=16
#SBATCH --job-name=snakemake
#SBATCH --partition campus-new
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=150

snakemake --cores $2 $1

```

After creating the shell script, and making it executable, run:

```
sbatch -o /path/to/log/file.log name_of_script.sh rule_name num_cores
```

## Log Files

For the `fimo` and `gene_x_tf rules`, separate log files are created and can be accessed in the `logs/` directory.

Other log files can be accessed in the `.snakemake/log` directory, although they are not very informative.

*Improved logging is in progress.*

## References

The custom MOTIF annotations are based off of an idea called *in silico ChIP-seq* introduced in [this paper](https://www.biorxiv.org/content/10.1101/2022.06.15.496239v1)
