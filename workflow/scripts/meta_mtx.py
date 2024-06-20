import scanpy as sc
import pandas as pd
import mudata as md
from scipy import io
import os
import json

def main():
    # Load data
    print('Loading data...')
    atac_sc_ad = md.read(snakemake.params.sc_atac)
    atac_meta_ad = md.read(snakemake.params.atac)
    seacell_label = snakemake.params.seacell_label
    seacells = atac_sc_ad.obs[seacell_label]
    group_variable = snakemake.params.cell_type_obs
    out_dir = snakemake.params.out_dir
    # find top cell type in each metacell
    top_ct = atac_sc_ad.obs[group_variable].groupby(
        seacells).value_counts().groupby(level=0, group_keys=False).head(1)
    atac_meta_ad.obs[group_variable] = top_ct[atac_meta_ad.obs_names].index.get_level_values(
        1)
    
    # Export data for edgeR
    print(f'Exporting data to {out_dir}...')
    
    # Create the directory, including all intermediate-level directories if necessary
    os.makedirs(str(out_dir), exist_ok = True)
    
    io.mmwrite(target=snakemake.output.out_x, a=atac_meta_ad.layers['raw'])
    pd.Series(atac_meta_ad.obs_names).to_csv(snakemake.output.out_obs, header = False)
    pd.Series(atac_meta_ad.var_names).to_csv(snakemake.output.out_var, header = False)
    atac_meta_ad.obs.to_csv(snakemake.output.out_metadata)
    
    
    with open (snakemake.output.json, 'w') as file:
        json.dump(snakemake.params.reference, file)
if __name__ == "__main__":
    main()
