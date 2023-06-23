if __name__ == "__main__":
    import argparse
    desc = "Selects primed- and lineage-specific ATAC-seq peaks based on user provided lineages and "
    
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--atac",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to ATAC single cell AnnData, common obs with RNA",
    )
    parser.add_argument(
        "--rna",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to RNA single cell AnnData, common obs with ATAC",
    )
    parser.add_argument(
        "--target",
        metavar="str",
        type=str,
        required=True,
        help="JSON array of target lineage and target cell type to convert to list",
    )
    parser.add_argument(
        "--start",
        type=str,
        required=True,
        help="Starting cell type name",
    )
    parser.add_argument(
        "--reference",
        metavar="str",
        type=str,
        required=True,
        help="JSON formatted string with reference lineages as keys and associated cell type as values to convert to dictionary",
    )
    parser.add_argument(
        "--min_logFC",
        type=float,
        default=-0.25,
        help="Float value of minimum logFC threshold for target lineage specific peaks",
    )
    parser.add_argument(
        "--max_logFC",
        type=float,
        default=0.25,
        help="Float value of max logFC threshold for reference lineages",
    )
    parser.add_argument(
        "--max_pval",
        type=float,
        default=0.1,
        help="Float value of max p-value threshold for gene-peak correlations",
    )
    parser.add_argument(
        "--min_corr",
        type=float,
        default=0.0,
        help="Float value of minimum correlation coefficient for gene-peak correlations",
    )

    args = parser.parse_args()


import scanpy as sc
import numpy as np
import json


def main(args):
    print(args)
    atac_ad = sc.read_h5ad(args.atac)
    rna_ad = sc.read_h5ad(args.rna)
    #set target lineage and cell type
    target = json.loads(args.target) 
    target_lineage = target[0]
    target_celltype = target[1]
    ref_lineages = json.loads(args.reference)
    start_celltype = args.start
    #set filtering parameters
    min_logFC = args.min_logFC
    max_pval = args.max_pval
    min_cor = args.min_cor
    max_logFC = args.max_logFC

    gene_peak_corrs = rna_ad.uns['gp_corrs']
        
    diff_peaks = dict()
    start_peaks = dict()
    for ref_lineage in ref_lineages.keys():
        ref_celltype = ref_lineages[ref_lineage]
        #differential accessibility compared to target cell type

    #select list of peaks which are positively correlated with gene expression and meet p value criteria
    correlated_peaks = pd.Series()
    for gene in gene_peak_corrs.keys():
        if type(gene_peak_corrs[gene]) == int:
            continue
        correlated_peaks = pd.concat([correlated_peaks, gene_peak_corrs[gene][np.logical_and(gene_peak_corrs[gene]['cor'] > min_cor,gene_peak_corrs[gene]['pval'] < max_pval)]])
    correlated_peaks = correlated_peaks.loc[:, ['cor', 'pval']]

    #convert open peak matrix to dataframe to access it by columns (celltype)
    open_peaks = pd.DataFrame(atac_ad.varm['OpenPeaks'].todense(), index = atac_ad.var_names, columns = atac_ad.uns['OpenPeaksCelltypes'])

    #filter differentially accessible peaks compared to target cell type to those which are open in the target celltype
    print(f"Removing peaks with logFC < {min_logFC}, correlation <= {min_corr} or pvalue >= ")
    for ref_lineage in diff_peaks.keys():
        #open peak filtering
        binary_filter_counts = open_peaks.T.loc[target_celltype, diff_peaks[ref_lineage].index].sum(axis=0)
        binary_filter_index = [True if ( x > 0) else False for x in binary_filter_counts]
        #filter out peaks below the minimum logFC cutoff
        diff_peaks[ref_lineage] = diff_peaks[ref_lineage][binary_filter_index].copy()
        diff_peaks[ref_lineage] = diff_peaks[ref_lineage][diff_peaks[ref_lineage]['logFC'] >= min_logFC]
        #subset to only peaks which are significantly correlated
        diff_peaks[ref_lineage] = diff_peaks[ref_lineage].loc[diff_peaks[ref_lineage].index.isin(correlated_peaks.index), :]

    #create set of unique peaksm which 
    unique_peaks = dict()
    unique_peaks[target_lineage] = pd.DataFrame()
    for ref_lineage in diff_peaks.keys():
        unique_peaks[target_lineage] = pd.concat([unique_peaks[target_lineage], diff_peaks[ref_lineage]]).sort_values("logFC", ascending = False).drop_duplicates('feature')

    #use Open peaks matrix to determine priming status - primed if open in HSCs, lineage specific if not
    unique_peaks[target_lineage]['primed'] = pd.Series([True if x == 1 else False for x in open_peaks.T.loc[start_celltype, unique_peaks[target_lineage].index].tolist()], unique_peaks[target_lineage].index)
    unique_peaks[target_lineage]['lineage_specific'] = pd.Series([ not x for x in unique_peaks[target_lineage]['primed']], index= unique_peaks[target_lineage].index)

    #Remove peaks which have positive fold changes above cutoff in any other lineage
    binary_filter_index = np.ones_like([1 if x < (max_logFC) else 0 for x in start_peaks[ref_1_lineage].loc[unique_peaks[target_lineage].index, 'logFC']])
    for ref_lineage in diff_peaks.keys():
        ref_list = [1 if x < (max_logFC) else 0 for x in start_peaks[ref_lineage].loc[unique_peaks[target_lineage].index, 'logFC']]
        binary_filter_index = np.bitwise_and(ref_list, binary_filter_index)

    #retain all peaks which are associated with a primed peak or lineage specific peak
    binary_filter_index = np.logical_or(np.logical_and(binary_filter_index, unique_peaks[target_lineage]['primed']),unique_peaks[target_lineage]['lineage_specific'])
    unique_peaks[target_lineage]= unique_peaks[target_lineage][binary_filter_index]
    print("Peak selection completed! Saving metadata ...")
    # add unique peaks to annData
    atac_ad.var[f'{target_lineage}_unique_peaks'] = pd.Series(atac_ad.var_names.isin(unique_peaks[target_lineage].index), index = atac_ad.var_names)
    # save unique peak results into .uns
    atac_ad.write(args.atac)
    print("Completed!")
    
    
if __name__ == "__main__":
    main(args)
