if __name__ == "__main__":
    import argparse
    desc = "Helper script to add differential accessibility output to anndata."

    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--atac",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to ATAC meta-cell AnnData, common obs with RNA",
    )
    parser.add_argument(
        "--start",
        type=str,
        required=True,
        help="string corresponding to cell type of the starting/stem cell type.",
    )
    parser.add_argument(
        "--reference",
        type=dict,
        required=True,
        help="dict of lineages matched to mature cell types.",
    )
    parser.add_argument(
        "--data_dir",
        type=str,
        required=True,
        help="directory with edgeR results.",
    )

    args = parser.parse_args()

import scanpy as sc
import muata as md
import pandas as pd
from itertools import product

def get_non_target_keys(my_dict, target):
    # Extract all keys
    all_keys = list(my_dict.keys())
    
    # Filter out the key specified by the user
    remaining_keys = [key for key in all_keys if key != target]
    
    # Return the remaining keys
    return remaining_keys
    
def main(args):
    # Arguments
    to_compare = args.to_compare.split('/')
    data_dir = args.data_dir

    # Load anndata
    atac_ad = md.read(args.atac)
    for key in args.reference:
        reference = get_non_target_keys(args.reference, key)
        to_compare = product(reference, [start, key])
        # Load diff results and save to anndata
        for ct_1, ct_2 in range(len(to_compare)):
    
            # Diff results
            df = pd.read_csv(data_dir + f"/{ct_1}_{ct_2}_diff_acc.tsv", sep="\t", index_col = "feature")
            #df.index = df['feature']
    
            # Add to anndata
            atac_ad.varm[f'{ct_1}_{ct_2}_diff_acc'] = df.loc[atac_ad.var_names, :]

    # Save anndata
    md.write(args.atac, atac_ad)


if __name__ == "__main__":
    main(args)
