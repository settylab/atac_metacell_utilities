if __name__ == "__main__":
    import argparse
    desc = "Compute correlations between RNA gene expression and ATAC accessibility"
    
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--atac",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to ATAC MC AnnData, common obs with RNA"
    )
    
    parser.add_argument(
        "--rna",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to RNA MC AnnData, common obs with RNA"
    )
    
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Path to output directory",
        default="results/gp_corrs/",
        metavar="directory"
    )
    args = parser.parse_args()

import pandas as pd
import numpy as np
import scanpy as sc
import SEACells
import pickle

def main(args):
    # Load data
    atac_ad = sc.read(args.atac)
    rna_ad = sc.read(args.rna)
    
    gene_set = rna_ad.var_names
    
    gp_corrs = SEACells.genescores.get_gene_peak_correlations(atac_int_ad, rna_int_ad,
                                                              path_to_gtf='../data/hg38.gtf',
                                                              span=100000, n_jobs=1,
                                                              gene_set=gene_set)
    
    peak_counts = SEACells.genescores.get_peak_counts(gp_corrs)
    
    # save files
    with open(args.outdir + 'gp_corr.pickle', 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    with open(args.outdir + 'peak_cts.pickle', 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

    
if __name__ == "__main__":
    main(args)