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
        "--to_compare",
        type=str,
        required=True,
        help="string of cell type comparisons for differential accessibility- groups of 2 are separated by semi-colons," +
        "comma separation within groups.",
    )
    parser.add_argument(
        "--data_dir",
        type=str,
        required=True,
        help="directory with edgeR results.",
    )

    args = parser.parse_args()

import scanpy as sc
import pandas as pd


def main(args):
    # Arguments
    to_compare = args.to_compare.split('/')
    data_dir = args.data_dir

    # Load anndata
    atac_ad = sc.read(args.atac)

    # Load diff results and save to anndata
    for i in range(len(to_compare)):
        # Cell types
        to_compare[i] = to_compare[i].split(',')

        # Diff results
        df = pd.read_csv(data_dir + f"/{to_compare[i][0]}_{to_compare[i][1]}_diff_acc.tsv", sep="\t")
        df.index = df['feature']

        # Add to anndata
        atac_ad.varm[f'{to_compare[i][0]}_{to_compare[i][1]}_diff_acc'] = df.loc[atac_ad.var_names, :]

    # Save anndata
    atac_ad.write(args.atac)


if __name__ == "__main__":
    main(args)
