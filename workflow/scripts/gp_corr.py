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
    
    parser.add_argument(
        "--gtf_file",
        type=str,
        help="Number of jobs for computing GP correlations",
        default="data/hg38.gtf",
        metavar="directory"
    )
    
    parser.add_argument(
        "--n_jobs",
        type=int,
        help="Number of jobs for computing GP correlations",
        default=1,
        metavar="int"
    )

    parser.add_argument(
        "--transcript_span",
        type=int,
        help="num of base pairs to add on longest transcript",
        default=100000,
        metavar="int"
    )
    
    parser.add_argument(
        "--test_set",
        help="Test on a subset of genes, n=--n_genes",
        action="store_true"
    )
    
    parser.add_argument(
        "--n_genes",
        type=int,
        help="num of base pairs to add on longest transcript",
        default=20,
        metavar="int"
    )
    
    parser.add_argument(
        "--min_corr",
        metavar="float",
        type=float,
        default=0.0,
        help="Minimum (excl) correlation value",
    )
    
    parser.add_argument(
        "--max_pval",
        metavar="float",
        type=float,
        default=0.1,
        help="Maximum (excl) p-value",
    )

    args = parser.parse_args()

import scanpy as sc
import SEACells
import pickle
import pandas as pd
from tqdm.auto import tqdm

def main(args):
    # Load data
    atac_ad = sc.read(args.atac)
    rna_ad = sc.read(args.rna)
    
    if args.test_set:
        gene_set = rna_ad.var_names[:args.n_genes]
    else:
        gene_set = rna_ad.var_names
    
    # Compute Gene-Peak Correlation scores
    gp_corrs = SEACells.genescores.get_gene_peak_correlations(atac_ad, rna_ad,
                                                              path_to_gtf=args.gtf_file,
                                                              span=args.transcript_span, n_jobs=args.n_jobs,
                                                              gene_set=gene_set)
    
    #peak_counts = SEACells.genescores.get_peak_counts(gp_corrs, max_pval=args.max_pval, min_corr=args.min_corr)
    peak_counts = SEACells.genescores.get_gene_peak_assocations(gp_corrs, pval_cutoff=args.max_pval, cor_cutoff=args.min_corr)
    # save files
    with open(args.outdir + '/gp_corr.pickle', 'wb') as handle:
        pickle.dump(gp_corrs, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    with open(args.outdir + '/peak_cts.pickle', 'wb') as handle:
        pickle.dump(peak_counts, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    gp_corr_df = pd.DataFrame()
    for gene in tqdm(gp_corr.keys()):
        if type(gp_corr[gene]) != int:
            gp_corr[gene]["gene"] = pd.Series(gene, gp_corr[gene].index)
            gp_corr_df = pd.concat([gp_corr_df, gp_corr[gene]])
    atac_ad.uns["gp_corrs"] = gp_corr_df
    
if __name__ == "__main__":
    main(args)
    