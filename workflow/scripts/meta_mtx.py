if __name__ == "__main__":
    import argparse
    desc = "Output ATAC metacell pseudo-bulk data as .mtx  and metadata as .csv for analysis in EdgeR."

    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--atac",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to ATAC MC AnnData, common obs with RNA",
    )
    parser.add_argument(
        "--sc_atac",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to ATAC single-cell AnnData, common obs with RNA",
    )
    parser.add_argument(
        "-o",
        metavar="path",
        type=str,
        required=True,
        help="Destination path for output files",
    )
    parser.add_argument(
        "--cell_type_obs",
        type=str,
        required=True,
        help="Cell type (or other column) name in .obs of single-cell anndata for grouping ",
    )
    parser.add_argument(
        "--seacell_label",
        type=str,
        default="SEACell",
        help="SEACell column name in .obs of single-cell anndata",
    )

    parser.add_argument(
        "--modality",
        type=str,
        default="atac",
        help="Change file labels when using rule with different modality.",
    )

    args = parser.parse_args()


import scanpy as sc
import pandas as pd
from scipy import io


def main(args):
    # Load data
    print('Loading data...')
    atac_sc_ad = sc.read(args.sc_atac)
    atac_meta_ad = sc.read(args.atac)
    seacell_label = args.seacell_label
    seacells = atac_sc_ad.obs[seacell_label]
    group_variable = args.cell_type_obs

    # find top cell type in each metacell
    top_ct = atac_sc_ad.obs[group_variable].groupby(
        seacells).value_counts().groupby(level=0, group_keys=False).head(1)
    atac_meta_ad.obs[group_variable] = top_ct[atac_meta_ad.obs_names].index.get_level_values(
        1)

    # Export data for edgeR
    print('Exporting data...')
    io.mmwrite(target=f"{args.o}/meta_{args.modality}_X.mtx", a=atac_meta_ad.layers['raw'])
    pd.Series(atac_meta_ad.obs_names).to_csv(f"{args.o}/meta_{args.modality}_cells.csv")
    pd.Series(atac_meta_ad.var_names).to_csv(f"{args.o}/meta_{args.modality}_peaks.csv")
    atac_meta_ad.obs.to_csv(f"{args.o}/meta_{args.modality}_metadata.csv")

if __name__ == "__main__":
    main(args)
