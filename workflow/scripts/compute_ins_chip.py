if __name__ == "__main__":
    import argparse
    desc = "Computes tf-peak correlations and in silico ChIP matrix"
    
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
        "--peak_tf",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to peak x TF AnnData",
    )
    
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Path to output directory",
        default="results/ins_chip/",
        metavar="directory"
    )
    
    parser.add_argument(
        "--verbose",
        help="Whether to print out info about each TF",
        action="store_true"
    
    )
    
    args = parser.parse_args()

    
import pandas as pd
import numpy as np
import scipy.io
from tqdm.auto import tqdm
from scipy.sparse import csr_matrix
import scanpy as sc
from anndata import AnnData
import numba as nb
import os

def filter_gene(ad, gene_list):
    has_expr = [g in gene_list for g in ad.var_names]
    
    frac = sum(has_expr) / len(has_expr)
    print(f"{frac:.2%} ({sum(has_expr)}) TFs have expresssion")
    
    ad = ad[:, has_expr]
    
    return ad

def mm_normalize(x):
    return (x - np.min(x)) / (np.max(x) - np.min(x))


@nb.njit
def mean1(a):
    n = len(a)
    b = np.empty(n)
    for i in range(n):
        b[i] = a[i].mean()
    return b


@nb.njit
def std1(a):
    n = len(a)
    b = np.empty(n)
    for i in range(n):
        b[i] = a[i].std()
    return b


@nb.njit
def c(a, b):
    m, k = b.shape

    mu_a = np.mean(a)
    mu_b = mean1(b)
    sig_a = np.std(a)
    sig_b = std1(b)
    
    out = np.empty((m, 1))

    for i in range(m):
        out[i] =  np.dot((a - mu_a),(b[i] - mu_b[i]))  / k / sig_a / sig_b[i]

    return out


def compute_corrs(tf_set, rna_mat, atac_mat, fimo_mat, verbose=True):
    corr_mat = pd.DataFrame(0.0, columns=tf_set, index=atac_mat.columns)
    for tf in tqdm(tf_set, total=len(tf_set)):
        relevant_peaks = fimo_mat.index[fimo_mat[tf] !=0]
        if verbose:
            print(f'{tf}: {len(relevant_peaks)} non-zero peaks')
        tf_expr = rna_mat[tf]
        atac_acc = atac_mat[relevant_peaks]

        # rank values
        tf_rank = scipy.stats.rankdata(tf_expr)
        peak_ranks = scipy.stats.rankdata(atac_acc, axis=0)  
        cor = c(tf_rank, peak_ranks.T)

        corr_mat.loc[relevant_peaks, tf] = cor[:,0]
    
    return corr_mat


def compute_ins_mat(tf_set, corr_mat, atac_mat, fimo_mat, verbose=True):
    final_out = pd.DataFrame(0.0,index=atac_mat.columns, columns=tf_set)
    for tf in tqdm(tf_set,total=len(tf_set)):
        pos_peaks = corr_mat.index[corr_mat[tf] > 0]

        if verbose:
            print(f'{tf}: {len(pos_peaks)} pos corr peaks')
    
        max_acc = atac_mat[pos_peaks].max(axis=0)
        
        corr_score = corr_mat.loc[pos_peaks,tf]
        fimo_score = fimo_mat.loc[pos_peaks, tf]
        fimo_score_n = fimo_score / fimo_score.max()

        mmn_fimo = mm_normalize(max_acc * fimo_score_n)
        pred = corr_score * mmn_fimo

        df = pd.DataFrame.from_dict({'corr': corr_score, 'max_accessibility': max_acc, 'fimo_score':fimo_score, 
                                     'fimo_score_n':fimo_score_n, 'score':pred})
        
        final_out.loc[pos_peaks, tf] = df.loc[pos_peaks, 'score']
        if verbose:
            print(f'completed for {tf}')

    return final_out
    
    
def main(args):
    # Load data
    atac_ad = sc.read(args.atac)
    rna_ad = sc.read(args.rna)
    peak_x_tf = sc.read(args.peak_tf)

    # make sure peaks are the same order:
    peak_x_tf = peak_x_tf[atac_ad.var_names, :].copy()
    
    # Filter genes for ones that have FIMO score
    gene_set = rna_ad.var_names
    peak_x_tf = filter_gene(peak_x_tf, gene_set)

    # grab matrices from anndatas
    rna_mat = rna_ad.to_df()
    atac_mat = atac_ad.to_df()
    fimo_mat = peak_x_tf.to_df()
    
    all_tfs = fimo_mat.columns
    
    # compute correlation and in silico chip matrix
    print('computing correlations...')
    corr_mat = compute_corrs(all_tfs,rna_mat, atac_mat, fimo_mat, verbose=args.verbose)
    print('Computing in silico chip score...')
    insc_mat = compute_ins_mat(all_tfs, corr_mat, atac_mat, fimo_mat, verbose=args.verbose)
    
    ## EXPORT DATA ##
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    
    # correlation matrix
    corr_outfile = args.outdir +'/tf_peak_corrs.mtx'
    print(f'writing tf-peak correlations to {corr_outfile}...')
    scipy.io.mmwrite(corr_outfile, csr_matrix(corr_mat.values), )    
    
    # in silico ChIP matrices
    ins_outfile = args.outdir +'/ins_chip.mtx'
    print(f'writing in silico ChIP matrix to {ins_outfile}...')
    scipy.io.mmwrite(ins_outfile, csr_matrix(insc_mat.values), )
    

    # TF and peak list
    print(f'saving metadata...')
    pd.Series(all_tfs).to_csv(args.outdir +'/tf_names.csv', index=False, header=['tf_name'])
    pd.Series(atac_mat.columns).to_csv(args.outdir +'/peak_names.csv', index=False, header=['peak_name'])
       

    print('Completed!')
    
    
if __name__ == "__main__":
    main(args)
