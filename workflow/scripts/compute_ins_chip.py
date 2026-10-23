import pandas as pd
import numpy as np
import scipy.io
from tqdm.auto import tqdm
from scipy.sparse import csr_matrix
import scanpy as sc
import numba as nb
import os
from itertools import product

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
        out[i] = np.dot((a - mu_a), (b[i] - mu_b[i])) / k / sig_a / sig_b[i]

    return out


def compute_corrs(tf_set, rna_mat, atac_mat, fimo_mat, verbose=True):
    corr_mat = pd.DataFrame(0.0, columns=tf_set, index=atac_mat.columns)
    for tf in tqdm(tf_set, total=len(tf_set)):
        relevant_peaks = fimo_mat.index[fimo_mat[tf] != 0]
        if verbose:
            print(f'{tf}: {len(relevant_peaks)} non-zero peaks')
        tf_expr = rna_mat[tf]
        atac_acc = atac_mat[relevant_peaks]

        # rank values
        tf_rank = scipy.stats.rankdata(tf_expr)
        peak_ranks = scipy.stats.rankdata(atac_acc, axis=0)
        cor = c(tf_rank, peak_ranks.T)

        corr_mat.loc[relevant_peaks, tf] = cor[:, 0]

    return corr_mat


def compute_ins_mat(tf_set, corr_mat, atac_mat, fimo_mat, verbose=True):
    final_out = pd.DataFrame(0.0, index=atac_mat.columns, columns=tf_set)
    for tf in tqdm(tf_set, total=len(tf_set)):
        pos_peaks = corr_mat.index[corr_mat[tf] > 0]

        if verbose:
            print(f'{tf}: {len(pos_peaks)} pos corr peaks')

        max_acc = atac_mat[pos_peaks].max(axis=0)

        corr_score = corr_mat.loc[pos_peaks, tf]
        fimo_score = fimo_mat.loc[pos_peaks, tf]
        fimo_score_n = fimo_score / fimo_score.max()

        mmn_fimo = mm_normalize(max_acc * fimo_score_n)
        pred = corr_score * mmn_fimo

        df = pd.DataFrame.from_dict({'corr': corr_score, 'max_accessibility': max_acc, 'fimo_score': fimo_score,
                                     'fimo_score_n': fimo_score_n, 'score': pred})

        final_out.loc[pos_peaks, tf] = df.loc[pos_peaks, 'score']
        if verbose:
            print(f'completed for {tf}')

    return final_out


def main():

    # Load data
    print('Load data...')
    atac_ad = sc.read(snakemake.config["anndata"]['atac'])
    rna_ad = sc.read(snakemake.config["anndata"]['rna'])
    atac_sc_ad = sc.read(snakemake.config["anndata"]['sc_atac'])
    rna_sc_ad = sc.read(snakemake.config["anndata"]['sc_rna'])
    # Reconstruction peak X TF
    peak_x_tf = sc.AnnData(atac_sc_ad.varm['FIMO'])
    peak_x_tf.var_names = atac_sc_ad.uns['FIMOColumns']
    peak_x_tf.obs_names = atac_sc_ad.var_names

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
    corr_mat = compute_corrs(all_tfs, rna_mat, atac_mat,
                             fimo_mat, verbose=snakemake.params['verbose'])
    print('Computing in silico chip score...')
    insc_mat = compute_ins_mat(
        all_tfs, corr_mat, atac_mat, fimo_mat, verbose=snakemake.params['verbose'])

    # Save results to ATAC anndata
    print('Saving results to anndata...')
    atac_sc_ad.varm['InSilicoChip_Corrs'] = csr_matrix(corr_mat.loc[atac_sc_ad.var_names, :].values)
    atac_sc_ad.varm['InSilicoChip'] = csr_matrix(insc_mat.loc[atac_sc_ad.var_names, :].values)
    atac_sc_ad.uns['InSilicoChipColumns'] = np.array(all_tfs)
    atac_sc_ad.write(snakemake.config['anndata']['sc_atac'])

    # Save targets to RNA anndata
    near_targets = dict()
    target_genes = rna_sc_ad[:,rna_sc_ad.var['upregulated_genes']].var_names
    near_target_genes = pd.DataFrame(0, index = rna_sc_ad.var_names, columns = target_genes )
    for chip_gene, gene in list(product(all_tfs, target_genes)):
        if chip_gene not in near_targets.keys():
            near_targets[chip_gene] = pd.DataFrame()
        near = insc_mat.loc[atac_ad.var['nearestGene'] == gene].copy()
        near.loc[:,'gene'] = pd.Series(data = gene, index =near.index)
        near = near.loc[near[chip_gene] >0,[chip_gene, 'gene']]
        near_targets[chip_gene] = pd.concat([near_targets[chip_gene], near])
    near_target_genes.loc[chip_gene, :]= [0 if x in near_targets[chip_gene]['gene'].unique() else 1 for x in near_target_genes.columns]
    print(near_target_genes.values)
    print(near_target_genes.shape)
    for col in rna_sc_ad.obs.select_dtypes(exclude=["number","bool","object_"]).columns:
        print(rna_sc_ad.obs[col].dtype)
        rna_sc_ad.obs[col] = rna_sc_ad.obs[col].astype(str).astype("category")
    rna_sc_ad.varm['TF_targets'] = near_target_genes.astype(int).values
    rna_sc_ad.uns['TF_targets_columns'] = list(target_genes.astype(str))
    rna_sc_ad.write_h5ad(snakemake.config['anndata']['sc_rna'])

    ## Create directory to signify completion
    os.makedirs(f"{str(snakemake.config['output'])}ins_chip/", exist_ok=True)
    with open(snakemake.output, 'a'):
        os.utime(snakemake.output, None)
    print('Completed!')


if __name__ == "__main__":
    main()
