if __name__ == "__main__":
    import argparse
    desc = "Compute open peaks and accessibility metrics -- the fraction of correlated open peaks"
    
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
        "-o",
        "--outdir",
        type=str,
        help="Path to output directory",
        default="results/accessibility/",
        metavar="directory"
    )
    
    parser.add_argument(
        "--gp_corr",
        metavar="directory",
        type=str,
        default="results/gp_corr/",
        help="Path to gene-peak correlation output directory",
    )
    
    parser.add_argument(
        "--embedding",
        type=str,
        help="`atac.obsm` field for nearest neighbor computation",
        default='X_svd',
        metavar="str"
    )
    
    parser.add_argument(
        "--op_pval",
        type=float,
        help="Nominal p-value cutoff for open peaks",
        default=1e-2,
        metavar="float"
    )
    parser.add_argument(
        "--read_len",
        type=int,
        help="Fragment length",
        default=147,
        metavar="int"
    )
    
    parser.add_argument(
        "--n_nbors",
        type=int,
        help="Number of neighbors for kNN",
        default=3,
        metavar="int"
    )
    
    # gene accessibility params
    parser.add_argument(
        "--gp_pval",
        type=float,
        help="Nominal p-value cutoff for gene-peak correlations",
        default=0.1,
        metavar="float"
    )
    
    parser.add_argument(
        "--min_corr",
        type=float,
        help="Correlation cuttoff",
        default=0.1,
        metavar="float"
    )
    
    parser.add_argument(
        "--min_peaks",
        type=int,
        help="Minimum number of peaks, inclusive",
        default=5,
        metavar="int"
    )
    
    args = parser.parse_args()

import scanpy as sc
import SEACells
import pickle


def main(args):
    # Load data
    atac_ad = sc.read(args.atac)
    
    with open(args.gp_corr + "/gp_corr.pickle", 'rb') as handle:
        gene_peak = pickle.load(handle)
        
    with open(args.gp_corr + "/peak_cts.pickle", 'rb') as handle:
        peak_counts = pickle.load(handle)
        
    # Determine open peaks
    SEACells.accessibility.determine_metacell_open_peaks(atac_ad, peak_set=None,low_dim_embedding=args.embedding,
                                                         pval_cutoff=args.op_pval, read_len=args.read_len,
                                                         n_neighbors=args.n_nbors)
    
    # compute gene accessibility
    high_reg_genes = peak_counts.index[peak_counts >= args.min_peaks]
    SEACells.accessibility.get_gene_accessibility(atac_ad, gene_peak_cors=gene_peak, gene_set=high_reg_genes,
                                                  pval_cutoff=args.gp_pval, cor_cutoff=args.min_corr)
    
    # grab 
    open_peaks = atac_ad.layers['OpenPeaks']
    accessibility = atac_ad.obsm['GeneAccessibility']
   
        
    # save files
    with open(args.outdir + '/open_peaks.pickle', 'wb') as handle:
        pickle.dump(open_peaks, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    with open(args.outdir + '/access_scores.pickle', 'wb') as handle:
        pickle.dump(accessibility, handle, protocol=pickle.HIGHEST_PROTOCOL)

    
if __name__ == "__main__":
    main(args)
