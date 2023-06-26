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
        help="Path to ATAC single-cell AnnData, common obs with RNA",
    )
    parser.add_argument(
        "--to_compare",
        type=str,
        required = True,
        help="string of cell type comparisons for differential accessibility- groups of 2 are separated by semi-colons, comma separation within groups.",
    )

    args = parser.parse_args()
import scanpy as sc
import pandas as pd

def main(args):
    to_compare = to_compare.split('/')
    atac_ad = sc.read(args.atac)
    for i in range(len(to_compare)):
        to_compare[i] = to_compare[i].split(',')
        df = pd.read_csv(data_dir + f"{to_compare[i][0]}_{to_compare[i][0]}_diff_acc.tsv", sep = "\t")
        atac_ad.varm[f'{to_compare[i][0]}_{to_compare[i][1]}_diff_acc'] = df.loc[atac_ad.var_names, :]
        atac_ad.uns[f'diff_acc_columns'] = df.columns
    atac_ad.write(args.atac)