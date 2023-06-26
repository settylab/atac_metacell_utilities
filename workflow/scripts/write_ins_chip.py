if __name__ == "__main__":
    import argparse
    desc = "add in-silico ChIP results to AnnDatas."
    
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
        "--datadir",
        type=str,
        required=True,
        help="Path to in silico ChIP results",
    )
    args = parser.parse_args()
    
import pandas as pd
import scanpy as sc
import numpy as np
from scipy.io import mmread
from itertools import product
def main(args):
    atac_ad = sc.read_h5ad(args.atac)
    rna_ad = sc.read_h5ad(args.rna)
    insc_mat = mmread(args.datadir+'/ins_chip.mtx')
    tf_names = pd.read_csv(args.datadir +'/tf_names.csv')
    tf_names = tf_names.iloc[:,0]
    insc_df = pd.DataFrame(insc_mat.todense(), index = atac_ad.var_names, columns = tf_names)
    print(insc_df)
    #print(f'writing in silico ChIP matrix to {ins_outfile}...')
    atac_ad.varm["in_silico_ChIP"]=insc_df.values
    atac_ad.uns['in_silico_ChIP_columns'] = list(insc_df.columns)
    print(f"writing in silico ChIP matrix to ATAC AnnData...")
    atac_ad.write_h5ad(args.atac)
    #save TF targets to 
    near_targets = dict()
    target_genes = rna_ad[:,rna_ad.var['upregulated_genes']].var_names
    near_target_genes = pd.DataFrame(0, index = rna_ad.var_names, columns = target_genes )
    for chip_gene, gene in list(product(tf_names, target_genes)):
        if chip_gene not in near_targets.keys():
            near_targets[chip_gene] = pd.DataFrame()
        near = insc_df.loc[atac_ad.var['nearestGene'] == gene].copy()
        near.loc[:,'gene'] = pd.Series(data = gene, index =near.index)
        near = near.loc[near[chip_gene] >0,[chip_gene, 'gene']]
        near_targets[chip_gene] = pd.concat([near_targets[chip_gene], near])
    near_target_genes.loc[chip_gene, :]= [0 if x in near_targets[chip_gene]['gene'].unique() else 1 for x in near_target_genes.columns]
    print(near_target_genes.values)
    print(near_target_genes.shape)
    for col in rna_ad.obs.select_dtypes(exclude=["number","bool","object_"]).columns:
        print(rna_ad.obs[col].dtype)
        rna_ad.obs[col] = rna_ad.obs[col].astype(str).astype("category")
    rna_ad.varm['TF_targets'] = near_target_genes.astype(int).values
    rna_ad.uns['TF_targets_columns'] = list(target_genes.astype(str))
    rna_ad.write_h5ad(args.rna)
    
    
if __name__ == "__main__":
    main(args)