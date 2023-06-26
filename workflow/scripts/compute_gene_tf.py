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
        "--rna",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to RNA MC AnnData, common obs with ATAC",
    )
    
    parser.add_argument(
        "--ins_chip",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to in-silico ChIP AnnData",
    )
    
    parser.add_argument(
        "--gp_corr",
        metavar="directory",
        type=str,
        default="results/gp_corr/",
        help="Path to gene-peak correlation output directory",
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
    
    parser.add_argument(
        "--min_peaks",
        metavar="int",
        type=int,
        default=2,
        help="Minimum num of significant peaks for gene to pass filtering, inclusive",
    )
    
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Path to output directory",
        default="results/gtf_mats/",
        metavar="directory"
    )
    
    args = parser.parse_args()

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import sparse
import pickle
from tqdm.auto import tqdm
from anndata import AnnData
    
def write_pickle(file_path, obj):    
    with open(file_path, 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)    

def filter_gene(ad, gene_list):
    has_expr = [g in gene_list for g in ad.var_names]
    
    frac = sum(has_expr) / len(has_expr)
    print(f"{frac:.2%} ({sum(has_expr)}) TFs have expresssion")
    
    ad = ad[:, has_expr]
    
    return ad

def count_genes_with_peaks(gene_list, gene_peak_scores):
    has_peak = [g in gene_peak_scores and isinstance(gene_peak_scores[g], pd.DataFrame)
                for g in gene_list]
    frac_peak = sum(has_peak) / len(has_peak)
    print(f"{frac_peak:.2%} ({sum(has_peak)}) genes have peak data")

def gene_tf_associations(gene_peak, peak_x_tf, 
                         min_corr=0.0, max_pval=0.1, min_peaks=2):
    
    tfs = peak_x_tf.var_names
    n_tfs=len(tfs)
    
    # initialize returns
    gene_tfs = {}
    n_gp_assoc = 0

    for gene, peak_df in tqdm(gene_peak.items(), 
                              desc="genes", 
                              total=len(gene_peak)):
        
        if gene not in gene_peak or not isinstance(peak_df, pd.DataFrame):
            # skip genes with no peaks
            continue

        sig_peaks = (peak_df['pval'] < max_pval) & (peak_df['cor'] > min_corr)

        if all(~sig_peaks):
            # no significant peaks
            continue
        
        # grab list of sig peaks correlated with this gene
        
        peaks = peak_df.loc[sig_peaks, :].index

        if len(peaks) >= min_peaks:

            # grab scores of all TFs for the sig peaks
            tf_association = peak_x_tf[peak_x_tf.obs_names.isin(peaks), :].X

            # convert scores to bool, 0 = False (peak x tf)
            tf_association = tf_association.toarray().astype(bool)

            # Grab TFs that are associated with a peak at least once 
            tf_idxs = np.where(np.any(tf_association, axis=0))[0]

            # tfs associated with each gene
            gene_tfs[gene] = {}

            for tf in tf_idxs:

                tf_name = peak_x_tf.var_names[tf]

                tf_mask = tf_association[:, tf] # subset the peak x tf scores to the relevant tfs
                gene_tfs[gene][tf_name] = peaks[tf_mask]

                # grab the peaks for each tf
                n_gp_assoc += 1
    
    # log number of genes that passed thresholds
    print(f"{len(gene_tfs):,} genes and {n_gp_assoc:,} gene-TF combinations with at least {min_peaks} peaks.")

    return gene_tfs


def compute_gene_tf_mat(atac_ad, peak_x_tf, gene_peak, gene_tfs):
    
    tfs = peak_x_tf.var_names

    # initiate three DFs
    gtf_sum = pd.DataFrame(0.0,index=gene_peak.keys(), columns=tfs)
    gtf_fimo, gtf_weighted = gtf_sum.copy(), gtf_sum.copy()
  
    atac_expr = pd.DataFrame(atac_ad.X.todense(), index=atac_ad.obs_names, 
                             columns=atac_ad.var_names)
    
    fimo_scores = pd.DataFrame(peak_x_tf.X.todense(), index=peak_x_tf.obs_names, 
                               columns=peak_x_tf.var_names)
    
    for gene, tf_dict in tqdm(gene_tfs.items(),total=len(gene_tfs)):
        for tf, peaks in tf_dict.items():

            # sum across cells for each peak
            sum_acc = atac_expr.loc[:, peaks].sum().sum()
            
            # Sum up FIMO scores for all peaks
            fimo_peaks = fimo_scores.loc[peaks, tf]
            fimo_sc = fimo_peaks.sum()
            
            # Compute weighted accessibility score
            weighted = atac_expr.loc[:, peaks].apply(lambda row: row.values * fimo_peaks.values, axis=1)
            
            gtf_sum.loc[gene, tf] = sum_acc
            gtf_fimo.loc[gene, tf] = fimo_sc
            gtf_weighted.loc[gene, tf] = weighted.sum().sum()
    
    return gtf_sum, gtf_fimo, gtf_weighted
            
def main(args):
    # Load data
    atac_ad = sc.read(args.atac)
    rna_ad = sc.read(args.rna)
    peak_x_tf_df = pd.DataFrame(atac_ad.varm["in_silico_ChIP"], index = atac_ad.var_names,  columns = atac_ad.uns['in_silico_ChIP_columns'])
    peak_x_tf = AnnData(peak_x_tf_df)
    peak_x_tf.X = sparse.csr_matrix(peak_x_tf.X) 
    
    gene_peak = dict(list(atac_ad.uns['gp_corrs'].groupby('gene')))
    
    # Filter genes for ones that have FIMO score
    gene_set = rna_ad.var_names
    peak_x_tf = filter_gene(peak_x_tf, gene_set)

    # Log genes with peaks
    count_genes_with_peaks(gene_set, gene_peak)
    
    # Create gene -> TF -> Peak dictionaries
    gene_tfs = gene_tf_associations(gene_peak, peak_x_tf, min_corr=args.min_corr,
                                    max_pval=args.max_pval, min_peaks=args.min_peaks)
    
    # Compute gene x TF matrices
    gtf_sum, gtf_fimo, gtf_weighted = compute_gene_tf_mat(atac_ad, peak_x_tf, gene_peak, gene_tfs)

    
    for mat in [('/sum', gtf_sum), ('/fimo', gtf_fimo), ('/weighted_sum', gtf_weighted)]:
        mat[1].to_csv(args.outdir + mat[0] + '.csv')

        # Pickle objects:
        write_pickle(args.outdir +  mat[0] + '.pickle', mat[1])
        
if __name__ == "__main__":
    main(args)
