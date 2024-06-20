import scanpy as sc
import mudata as md
import pandas as pd
import numpy as np
from tqdm.auto import tqdm


def main():

    # Load data
    print('Loading anndata...')
    atac_ad = md.read(snakemake.params["sc_atac"])
    atac_meta_ad = md.read(snakemake.params["meta_atac"])

    # Start cell type
    start_celltype = snakemake.params["start"]
    # Reference lineages
    reference = snakemake.params["reference"]

    # set filtering parameters
    min_logFC = snakemake.params["min_logFC"]
    max_logFC = snakemake.params["max_logFC"]
    max_pval = snakemake.params["max_pval"]
    min_corr = snakemake.params["min_corr"]


    print("Loading gene-peak corrs...")
    gene_peak_corrs = dict(list(atac_meta_ad.uns['gp_corrs'].groupby('gene')))
    for target_lineage in reference:
        ref_lineages = reference.copy()
        target_celltype = ref_lineages.pop[target_lineage]
        
        print("Loading differential results")
        diff_peaks = dict()
        start_peaks = dict()
        for ref_lineage in ref_lineages.keys():
            ref_celltype = ref_lineages[ref_lineage]
    
            # Differential accessibility compared to target cell type
            diff_peaks[ref_lineage] = atac_meta_ad.varm[f'{ref_celltype}_{target_celltype}_diff_acc']
            # diff_peaks[ref_lineage]['feature'] = diff_peaks[ref_lineage].index
    
            # Differential accessibility compared to stem cell type
            start_peaks[ref_lineage] = atac_meta_ad.varm[f'{ref_celltype}_{start_celltype}_diff_acc']
            # start_peaks[ref_lineage]['feature'] = diff_peaks[ref_lineage].index
    
        # Select list of peaks which are positively correlated with gene expression and meet p value criteria
        print("Removing peaks not positively correlated with gene expression...")
        correlated_peaks = pd.Series(dtype=object)
        for gene in tqdm(gene_peak_corrs.keys()):
            if isinstance(gene_peak_corrs[gene], int):
                continue
            res = gene_peak_corrs[gene]
            iter_peaks = res.index[(res['cor'] >= min_corr) & (res['pval'] < max_pval)]
            correlated_peaks = iter_peaks.union(correlated_peaks)
        filtered_peaks = correlated_peaks.unique()
        print(f'Retained {len(filtered_peaks)} peaks')
    
    
        # Peaks that are open in the target lineage
        print('Retain peaks that are open in the target lineage')
        open_peaks = atac_ad.varm['OpenPeaks']
        filtered_peaks = filtered_peaks[open_peaks.loc[filtered_peaks, [target_celltype]].sum(axis=1) > 0]
        print(f'Retained {len(filtered_peaks)} peaks')
    
    
        print(f"Retain peaks with greater accessibility in {target_lineage} compared to other lineages")
        retain_peaks = pd.Series(False, index=filtered_peaks)
        for ref in diff_peaks.keys():
            diff = diff_peaks[ref].loc[filtered_peaks]
            retain = filtered_peaks[(diff['logFC'] > min_logFC) & ~np.isnan(diff['groupA_N'])]
            retain_peaks[retain] = True
        filtered_peaks = filtered_peaks[retain_peaks]
        print(f'Retained {len(filtered_peaks)} peaks')
    
    
        print("Remove peaks with increased accessibility in reference cell types compared to stem cells")
        retain_peaks = pd.Series(False, index=filtered_peaks)
        for ref in start_peaks.keys():
            diff = start_peaks[ref].loc[filtered_peaks]
            retain = filtered_peaks[(diff['logFC'] < max_logFC)]
            retain_peaks[retain] = True
        filtered_peaks = filtered_peaks[retain_peaks]
        print(f'Retained {len(filtered_peaks)} peaks')
    
    
        # Primed and lineage specific peaks
        primed_peaks = filtered_peaks[open_peaks.loc[filtered_peaks, start_celltype] > 0]
        lin_spec_peaks = filtered_peaks[open_peaks.loc[filtered_peaks, start_celltype] == 0]
    
        print("Peak selection completed! Adding results to anndata ...")
        atac_ad.var[f'{target_lineage}_primed'] = atac_ad.var_names.isin(primed_peaks)
        atac_ad.var[f'{target_lineage}_lineage_specific'] = atac_ad.var_names.isin(lin_spec_peaks)

    print('Saving')
    md.write(snakemake.input["sc_atac"], atac_ad.write)

    # CReate directory to mark completion
    import os
    os.makedirs(str(snakemake.output["out_dir"]), exist_ok=True)

    print("Completed!")

if __name__ == "__main__":
    main()
