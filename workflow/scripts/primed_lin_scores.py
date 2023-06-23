#Perform actual imputation only on the peaks associated with genes of interest. I don't think it needs to be a dataframe as long as the shape is correct
if __name__ == "__main__":
    import argparse
    import pickle
    desc = "Computes weighted averages of imputed ATAC-seq accessibility for primed and lineage specific peaks."
    
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
        "--max_pval",
        type=float,
        default=0.1
        help="Float value of max p-value threshold for gene-peak correlations",
    )
    parser.add_argument(
        "--min_corr",
        type=float,
        default=0.0
        help="Float value of minimum correlation coefficient for gene-peak correlations",
    )

    args = parser.parse_args()

import scanpy as sc
import numpy as np
from tqdm.auto import tqdm
from scipy import sparse

def filter_gp_corrs(gp_corrs, gene_list, peak_df, min_cor, max_pval):
    filtered_gp_corrs = dict()
    for gene in list(gp_corrs.keys()):
        if type(gp_corrs[gene]) == int:
            continue
        elif gene not in gene_list:
            continue
        else:
            #restrict to only peaks in unique peak set
            filtered_gp_corrs[gene] = gp_corrs[gene][gp_corrs[gene].index.isin(peak_df.index)].copy()
            #restrict to peaks which pass correlation filter
            filtered_gp_corrs[gene] = filtered_gp_corrs[gene][filtered_gp_corrs[gene]['cor']>min_cor]
            #restrict to peaks which pass pvalue filter
            filtered_gp_corrs[gene] = filtered_gp_corrs[gene][filtered_gp_corrs[gene]['pval']<0.1]
            category  = pd.Series("other", index = filtered_gp_corrs[gene].index)
            category[np.where(peak_df.loc[category.index, 'primed'] == True)[0]] = 'primed'
            category[np.where(peak_df.loc[category.index,'lineage_specific'] == True)[0]] = 'lineage_specific'
            filtered_gp_corrs[gene]['category'] = category
    return filtered_gp_corrs

def _dot_func(x, y):
    return x.dot(y)

def impute_data_with_weights(ad, data):
    res = _dot_func(ad.obsp['ImputeWeights'], data)
    return pd.DataFrame(res, index=ad.obs_names, columns=ad.var_names)

def _concat_imputed_peak_values(ad, peak_df, gp_corrs):
    var_index = gp_corrs[gene][gp_corrs[gene]['category'] == key].index
    df = impute_data_with_weights(ad[:, var_index], pd.DataFrame(ad[:, var_index].layers['tf_idf'].todense(), columns = var_index, index = ad.obs_names))
    peak_df = pd.concat([peak_df, df] , axis = 1)

def _concat_weighted_peak_avg(gene_df, peak_df, gp_corrs, gene ):
    lineage_index = df.index
    filtered_peak_index = gp_corrs[gene][gp_corrs[gene]['category'] == key].index
    if len(filtered_peak_index) == 0:
        gene_df = pd.concat([gene_df, pd.Series(0.0, index = lineage_index, name = gene)], axis = 1)
    elif len(filtered_peak_index) == 1:
        gene_df = pd.concat([gene_df, pd.Series(peak_df.loc[lineage_index, filtered_peak_index].to_numpy().ravel(), index = lineage_index, name = gene)], axis = 1)
    else:
        gene_df = pd.concat([gene_df, pd.Series(np.average(peak_df.loc[lineage_index, filtered_peak_index], axis = 1, weights = gp_corrs[gene].loc[filtered_peak_index,'cor'].to_numpy()), index = lineage_index, name = gene)], axis = 1)


def main(args):
    atac_ad = sc.read_h5ad(args.atac)
    rna_ad = sc.read_h5ad(args.rna)

    gene_peak_corrs = rna_ad.uns['gp_corrs']

    target_genes = rna_ad[:, rna_ad.var['upregulated_genes']].index
    control_genes = rna_ad[:, rna_ad.var['control_genes']].index
    unique_peaks = atac_ad.var.loc[rna_ad.var[f'{target_lineage}_unique_peaks'], ['primed', 'lineage_specific']]

    atac_ad.uns['filtered_gp_corrs'] = filter_gp_corrs(gene_peak_corrs, target_genes, unique_peaks, args.min_cor, args.max_pval)
    atac_ad.uns['control_gp_corrs'] = filter_gp_corrs(gene_peak_corrs, control_genes, unique_peaks, args.min_cor, args.max_pval)

    imputed_peak_matrix = dict()
    for key in ['primed', 'lineage_specific']:
        imputed_peak_matrix[key] = pd.DataFrame(index = rna_ad.obs_names)
        #get imputed values for peaks associated with target genes
        for gene in tqdm(rna_ad.uns['filtered_gp_corrs'].keys()):
            _concat_imputed_peak_values(atac_ad, imputed_peak_matrix[key], rna_ad.uns['filtered_gp_corrs'])
        #get imputed values for peaks associated with control genes
        for gene in tqdm(rna_ad.uns['control_gp_corrs'].keys()):
            _concat_imputed_peak_values(atac_ad, imputed_peak_matrix[key], rna_ad.uns['control_gp_corrs'])
        #remove duplicates
        imputed_peak_matrix[key] = imputed_peak_matrix[key].T.drop_duplicates().T
    atac_ad.var['primed'] = pd.Series(0, index = atac_ad.var_names)
    atac_ad.var['lineage_specific'] = pd.Series(0, index = atac_ad.var_names)
    for gene in rna_ad.uns['filtered_gp_corrs']:
        index = rna_ad.uns['filtered_gp_corrs'][rna_ad.uns['filtered_gp_corrs']['category'] == "primed"].index
        atac_ad.var.loc[index, "primed"] = 1
        index = rna_ad.uns['filtered_gp_corrs'][rna_ad.uns['filtered_gp_corrs']['category'] == "lineage_specifc"].index
        atac_ad.var.loc[index, "lineage_specific"] = 1
        
    # weighted average
    atac_gene_dfs = dict()
    lineage_name = target_lineage+"_lineage"
    lineage_mask = rna_ad.obs[lineage_name].astype(bool)
    warnings.filterwarnings(category=FutureWarning, action = 'ignore')
    for key in imputed_peak_matrix.keys():
        lineage_index =rna_ad.obs_names
        atac_gene_dfs[key] = pd.DataFrame(index = lineage_index) 
        for gene in rna_ad.uns['filtered_gp_corrs'].keys():
            _concat_weighted_peak_avg(atac_gene_dfs[key], imputed_peak_matrix[key], rna_ad.uns['filtered_gp_corrs'], gene)
        for gene in rna_ad.uns['control_gp_corrs'].keys():
            _concat_weighted_peak_avg(atac_gene_dfs[key], imputed_peak_matrix[key], rna_ad.uns['control_gp_corrs'], gene)
        atac_gene_dfs[key] = atac_gene_dfs[key].fillna('0', axis = 1)  
    for key in atac_gene_dfs.key():
        gene_scores = pd.DataFrame(0.0, index = rna_ad.obs_names, columns = rna_ad.var_names)
        gene_scores.loc[:,atac_gene_dfs[key].columns] = atac_gene_dfs[key].values
        rna_ad.layers[f'{key}_scores'] = sparse.csr_matrix(gene_scores.values)
    rna_ad.write(args.rna)
    atac_ad.write(args.atac)
    
if __name__ == "__main__":
    main(args)