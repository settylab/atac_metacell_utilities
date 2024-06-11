if __name__ == "__main__":
    import argparse
    desc = "Adds Chromvar results to ATAC annData."

    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--sc_rna",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to RNA single-cell AnnData. Can also be path the the RNA modality of a MuData object.",
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
import mudata as md

def main(args):
    print('Loading and then saving chromVAR results to anndata..')
    rna_ad = md.read(args.sc_rna)
    deviations = pd.read_csv(args.input + "/deviations.csv", index_col=0).T
    deviations.index = rna_ad.obs_names
    rna_ad.obsm["chromVAR_deviations"] = deviations
    md.write(args.sc_rna, rna_ad)


if __name__ == "__main__":
    main(args)
