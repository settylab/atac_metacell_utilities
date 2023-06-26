if __name__ == "__main__":
    import argparse
    desc = "Adds Chromvar results to ATAC annData."
    
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--atac",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to ATAC single-cell AnnData",
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to ChromVAR results for input.",
    )
    

    args = parser.parse_args()

import scanpy as sc
import pandas as pd

def main(args)
    atac_ad = sc.read_h5ad(args.sc_atac)
    deviations_csv = pd.read_csv(args.input + "/deviations.csv", index_col = 0)
    atac_ad.obsm["chromVAR_deviations"] = deviations_csv.values
    atac_ad.uns["chromVAR_deviations_columns"] = deviations_csv.columns
