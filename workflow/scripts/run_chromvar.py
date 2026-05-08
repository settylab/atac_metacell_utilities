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
        help="Path to RNA single-cell AnnData",
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


def main(args):
    print('Loading and then saving chromVAR results to anndata..')
    rna_ad = sc.read_h5ad(args.sc_rna)
    deviations = pd.read_csv(args.input + "/deviations.csv", index_col=0).T
    
    # Match cells between ATAC (deviations) and RNA (rna_ad)
    # Find common cells
    common_cells = deviations.index.intersection(rna_ad.obs_names)
    
    if len(common_cells) == 0:
        raise ValueError(
            f"No common cells found between chromVAR deviations ({len(deviations)} cells) "
            f"and RNA AnnData ({len(rna_ad)} cells). "
            f"Deviations index sample: {deviations.index[:5].tolist()}, "
            f"RNA obs_names sample: {rna_ad.obs_names[:5].tolist()}"
        )
    
    print(f'Found {len(common_cells)} common cells out of {len(deviations)} ATAC cells and {len(rna_ad)} RNA cells')
    
    if len(common_cells) < len(deviations):
        print(f'Warning: {len(deviations) - len(common_cells)} ATAC cells not found in RNA data')
    
    if len(common_cells) < len(rna_ad):
        print(f'Warning: {len(rna_ad) - len(common_cells)} RNA cells not found in ATAC data')
    
    # Subset deviations to common cells and align with RNA AnnData order
    deviations_aligned = pd.DataFrame(
        index=rna_ad.obs_names,
        columns=deviations.columns,
        dtype=float
    )
    deviations_aligned.loc[common_cells] = deviations.loc[common_cells]
    
    # Fill missing cells with NaN (or zeros if preferred)
    # deviations_aligned = deviations_aligned.fillna(0)  # Uncomment if you want zeros instead of NaN
    
    rna_ad.obsm["chromVAR_deviations"] = deviations_aligned
    rna_ad.write(args.sc_rna)
    print(f'Successfully saved chromVAR deviations to {args.sc_rna}')


if __name__ == "__main__":
    main(args)
