import scanpy as sc
import numpy as np
import pandas as pd
from tqdm.auto import tqdm


def _dot_func(x, y):
    return x.dot(y)

def impute_data_with_weights(ad, data):
    res = _dot_func(ad.obsp['ImputeWeights'], data)
    return pd.DataFrame(res, index=ad.obs_names, columns=ad.var_names)


def main():

    # Load params
    print('Load data...')
    atac_sc_ad = sc.read_h5ad(snakemake.input["sc_atac"])
    atac_meta_ad = sc.read_h5ad(snakemake.input["meta_atac"])
    rna_sc_ad = sc.read_h5ad(snakemake.input["sc_rna"])
    target_lineage = snakemake.params["target"][0]
    min_corr = snakemake.params["min_corr"]
    max_pval = snakemake.params["max_pval"]
    min_peaks = snakemake.params["min_peaks"]

    # Reconstruct gene peak correlations
    gene_peak_corrs = dict(list(atac_meta_ad.uns['gp_corrs'].groupby('gene')))

    # Peak accessibility imputation
    print('Peak accessibility imputation...')
    lin_peaks = atac_sc_ad.var_names[
        (atac_sc_ad.var[f'{target_lineage}_primed']) |
        (atac_sc_ad.var[f'{target_lineage}_lineage_specific'])
    ]
    imputed_peak_matrix = impute_data_with_weights(atac_sc_ad[:, lin_peaks],
                                                   pd.DataFrame(atac_sc_ad[:, lin_peaks].layers['tf_idf'].todense(),
                                                                columns=lin_peaks, index=atac_sc_ad.obs_names))
    atac_sc_ad.var[f'{target_lineage}_all'] = pd.Series(atac_sc_ad.var_names.isin(lin_peaks), index = atac_sc_ad.var_names, dtype = bool)
    # Scores
    from scipy.sparse import csr_matrix
    acc_scores = dict()
    for key in ['primed', 'lineage_specific', 'all']:
        print(f'Computing {key} scores')

        # Imputead accessibility
        lin_peaks = atac_sc_ad.var_names[atac_sc_ad.var[f'{target_lineage}_{key}']]
        imp_peaks = imputed_peak_matrix[lin_peaks]

        # Scores per gene
        acc_scores[key] = pd.DataFrame(
            0.0, index=atac_sc_ad.obs_names, columns=rna_sc_ad.var_names)
        for gene in tqdm(rna_sc_ad.var_names):
            if gene not in gene_peak_corrs.keys():
                continue

            # Gene peak correlations
            res = gene_peak_corrs[gene]
            if isinstance(res, int):
                continue

            gene_peaks = res.index[(res['cor'] > min_corr)
                                   & (res['pval'] < max_pval)]
            if len(gene_peaks) < min_peaks:
                continue

            key_peaks = gene_peaks[gene_peaks.isin(imp_peaks.columns)]
            if len(key_peaks) == 0:
                continue

            acc_scores[key].loc[:, gene] = imp_peaks[key_peaks].dot(res.loc[key_peaks, 'cor']) \
                / np.sum(res.loc[key_peaks, 'cor'])

        # Save results to rna andata
        rna_sc_ad.layers[f'{target_lineage}_{key}'] = csr_matrix(acc_scores[key])
 
    # Save results
    print('Saving anndata')
    rna_sc_ad.write_h5ad(snakemake.input["sc_rna"])

    # CReate directory to mark completion
    import os
    os.makedirs(str(snakemake.output["out_dir"]), exist_ok=True)


if __name__ == "__main__":
    main()
