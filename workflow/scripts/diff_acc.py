import scanpy as sc
import mudata as md
import pandas as pd
from itertools import product
from pathlib import Path

def get_non_target_celltypes(my_dict, target):
    # Extract all keys
    all_keys = list(my_dict.keys())
    
    # Filter out the key specified by the user
    remaining_celltypes = [my_dict[key] for key in all_keys if key != target]
    
    # Return the remaining keys
    return remaining_celltypes
    
def main():
    # Arguments
    data_dir = snakemake.input.out_dir
    reference = snakemake.params.reference
    start = snakemake.params.start
    # Load anndata
    atac_ad = md.read(snakemake.params.atac)
    for key in reference:
        ref = get_non_target_celltypes(reference, key)
        to_compare = product(ref, [start, reference[key]])
        # Load diff results and save to anndata
        for ct_1, ct_2 in to_compare:
    
            # Diff results
            df = pd.read_csv(data_dir + f"/{ct_1}_{ct_2}_diff_acc.tsv", sep="\t", index_col = "feature")
            #df.index = df['feature']
    
            # Add to anndata
            atac_ad.varm[f'{ct_1}_{ct_2}_diff_acc'] = df.loc[atac_ad.var_names, :]

    # Save anndata
    md.write(snakemake.params.atac, atac_ad)
    Path(f'{data_dir}/.diff_acc_update').touch()
if __name__ == "__main__":
    main()
