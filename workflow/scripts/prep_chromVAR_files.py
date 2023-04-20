if __name__ == "__main__":
    import argparse
    desc = "Filters and binarizes chip matrix for chromvar"
    
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--sc_atac",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to ATAC single-cell AnnData",
    )
    
    parser.add_argument(
        "--ins_chip_dir",
        metavar="directory",
        type=str,
        required=True,
        help="Path to in silico ChIP output directory",
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
from scipy.sparse import csr_matrix, issparse  
import scanpy as sc
    
def filt_and_binarize(ins_chip_mat, min_chip_score, min_peak_hits):
    print(f'Filtering for TFs with insilico chip score > {min_chip_score} in at least {min_peak_hits} peaks')
    ins_chip_mat[ins_chip_mat <= min_chip_score] = 0.0
    ins_chip_bin = (ins_chip_mat != 0.0).astype(int)
    filter_bin = ins_chip_bin.loc[:,ins_chip_bin.sum() >= min_peak_hits]

    return filter_bin


def main(args):
    insc_mat = scipy.io.mmread(args.ins_chip_dir + '/ins_chip.mtx')
    tf_names = pd.read_csv(args.ins_chip_dir + '/tf_names.csv')['tf_name']
    peak_names = pd.read_csv(args.ins_chip_dir + '/peak_names.csv')['peak_name']
    
    insc_df = pd.DataFrame(insc_mat.A, columns=tf_names, index=peak_names)
    
    binary_mat = filt_and_binarize(insc_df, args.min_chip, args.min_peak_hits)
    print(f'Filtered {len(tf_names)} down to {binary_mat.shape[1]} TFs for chromVAR')

    # Binary for chromvar
    
    print(f'writing binary in silico ChIP matrix and TF names...')
    #export binary matrix for chromVAR
    scipy.io.mmwrite(args.outdir +'/bin_ins_chip.mtx', csr_matrix(binary_mat.values))
    # export tf names
    pd.Series(binary_mat.columns).to_csv(args.outdir +'/chromvar_tf_names.csv', index=False, header=['tf_name'])
    
    
    sc_atac_ad = sc.read(args.sc_atac)
    sc_atac_ad = sc_atac_ad[:,peak_names]
    
    if ~issparse(sc_atac_ad.X):
        sc_atac_ad.X = csr_matrix(sc_atac_ad.X)
    
    scipy.io.mmwrite(args.outdir +'/sc_atac_counts.mtx', sc_atac_ad.X)
    
    if 'nFrags' not in sc_atac_ad.obs.columns:
        raise KeyError("'nFrags' not in scATAC .obs")    
    sc_atac_ad.obs['nFrags'].to_csv(args.outdir + '/sc_nfrags.csv', index=True, header=['nFrags'])                             
                                 
    
if __name__ == "__main__":
    main(args)