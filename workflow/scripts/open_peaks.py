import scanpy as sc
import mudata as md
import SEACells
import pandas as pd
import numpy as np


def main():

    # Load data
    print('Loading data...')
    atac_sc_ad = md.read(snakemake.params["sc_atac"])
    group_variable = snakemake.params['cell_type_obs']

    # Open peaks per cell type
    # Aggregate counts per cell type
    print('Aggregating counts per cell type')
    ct_ad = SEACells.core.summarize_by_SEACell(atac_sc_ad,
                                               group_variable, summarize_layer='X')
    ct_ad.obs['n_counts'] = np.ravel(ct_ad.X.sum(axis=1))
    # SVD for filler - only used for identifying itself as nearest neighbor
    sc_svd = pd.DataFrame(atac_sc_ad.obsm['X_svd'], index=atac_sc_ad.obs_names)
    ct_ad.obsm['X_svd'] = sc_svd.groupby(
        atac_sc_ad.obs[group_variable]).mean().loc[ct_ad.obs_names, :]

    # Open peaks
    # Determine open peaks in each metacell
    print('Open peaks in each cell type')
    SEACells.accessibility.determine_metacell_open_peaks(
        ct_ad, n_neighbors=1, n_jobs=1)

    # Update anndata
    print('Update and saving anndata')
    atac_sc_ad.varm['OpenPeaks'] = pd.DataFrame(ct_ad.layers['OpenPeaks'].T,
                                                columns=ct_ad.obs_names, index=atac_sc_ad.var_names)

    # Save
    md.write(snakemake.input["sc_atac"],  atac_sc_ad)

    # CReate directory to mark completion
    import os
    os.makedirs(str(snakemake.output["out_dir"]), exist_ok=True)


if __name__ == "__main__":
    main()
