if __name__ == "__main__":
    import argparse
    desc = "Creates three versions of gene x TF matrices"
    
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
        "--ins_chip",
        metavar="sparse matrix, MTX file",
        type=str,
        required=True,
        help="Path to in silico ChIP matrix, peak x TF",
    )
    
    parser.add_argument(
        "--tf_names",
        metavar="TF names",
        type=str,
        required=True,
        help="Path to in silico ChIP matrix colnames",
    )
    
    parser.add_argument(
        "--peak_names",
        metavar="peak name list",
        type=str,
        required=True,
        help="Path to peak list",
    )
    
    parser.add_argument(
        "--min_chip",
        type=float,
        help="minimum in silico ChIP score to pass filtering",
        default=0.15,
        metavar="float"
    )
    
    parser.add_argument(
        "--min_peak_hits",
        type=int,
        help="minimum TFs passing in silico ChIP filters",
        default=30,
        metavar="int"
    )
    
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Path to output directory",
        default="results/chromvar/",
        metavar="directory"
    )


    args = parser.parse_args()

import pandas as pd
import numpy as np
import scipy.io
from scipy.sparse import csr_matrix    
    
def filt_and_binarize(ins_chip_mat, min_chip_score, min_peak_hits):
    print(f'Filtering for TFs with insilico chip score > {min_chip_score} in at least {min_peak_hits} peaks')
    ins_chip_mat[ins_chip_mat <= min_chip_score] = 0.0
    ins_chip_bin = (ins_chip_mat != 0.0).astype(int)
    filter_bin = ins_chip_bin.loc[:,ins_chip_bin.sum() >= min_peak_hits]

    return filter_bin


def main(args):
    insc_mat = scipy.io.mmread(args.ins_chip)
    insc_df = pd.DataFrame(insc_mat.A, columns=args.tf_names, index=args.peak_names)
    sc_atac_ad = args.sc_atac
    
    binary_mat = filt_and_binarize(insc_df, args.min_chip, args.min_peak_hits)
    print(f'Filtered {len(args.tf_names)} down to {binary_mat.shape[1]} TFs for chromVAR')

    # Binary for chromvar
    bin_ins_outfile = args.outdir +'/bin_ins_chip.mtx'
    if args.verbose:
        print(f'writing binary in silico ChIP matrix and TF names...')
    #export binary matrix for chromVAR
    scipy.io.mmwrite(bin_ins_outfile, csr_matrix(binary_mat.values))
    # export tf names
    pd.Series(binary_mat.columns).to_csv(args.outdir +'/chromvar_tf_names.csv', index=False, header=['tf_name'])
    
    if 'nFrags' not in sc_atac_ad.obs.columns:
        raise KeyError("'nFrags' not in scATAC .obs")    
    sc_atac_ad.obs['nFrags'].to_csv(args.outdir + '/sc_nfrags.csv', index=True, header=['nFrags']                             
                                 
    
if __name__ == "__main__":
    main(args)