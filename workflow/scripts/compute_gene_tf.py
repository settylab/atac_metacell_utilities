import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix

def main():
    # Load data
    print('Loading data')
    atac_sc_ad = sc.read(snakemake.config['anndata']['sc_atac'])
    atac_mc_ad = sc.read(snakemake.config['anndata']['atac'])
    rna_sc_ad = sc.read(snakemake.config['anndata']['sc_rna'])

    # Assemble necessary matrices
    print('Computing gene TF relationships')

    # Gene X Peak matrix
    gp = atac_mc_ad.uns['gp_corrs'].copy()
    gp = gp.loc[(gp['cor'] > snakemake.params['min_corr']) & (gp['pval'] < snakemake.params['max_pval']), :]
    peak_counts = gp.groupby('gene').apply(len)
    use_genes = peak_counts[peak_counts >= snakemake.params['min_peaks']].index
    gp.loc[gp['gene'].isin(use_genes), :]

    # Convert to categoricals
    gp['peaks'] = gp.index.values
    gp['peaks'] = pd.Categorical(gp['peaks'], atac_sc_ad.var_names)
    gp['gene'] = pd.Categorical(gp['gene'], rna_sc_ad.var_names)
    
    # Binarize gene peak correlations
    gp_binary = gp.copy()
    gp_binary['cor'] = (gp_binary['cor'] != 0).astype(int)
    
    # Construct sparse matrix
    data = gp_binary['cor'].values
    row_ind = gp_binary['gene'].values.codes
    col_ind = gp_binary['peaks'].values.codes
    M = len(gp_binary['gene'].values.categories)
    N = len(gp_binary['peaks'].values.categories)
    gene_X_peak = csr_matrix((data, (row_ind, col_ind)), [M, N])

    # Peak X TF matrix
    peak_X_tf = atac_sc_ad.varm['InSilicoChip']

    # Gene X TF matrix
    gene_X_tf = gene_X_peak.dot(peak_X_tf)
    
    rna_sc_ad.varm['geneXTF'] = pd.DataFrame(gene_X_tf.todense(),
        index=rna_sc_ad.var_names,
        columns=atac_sc_ad.uns['InSilicoChipColumns'])
    
    # Saving resuls
    print('Saving results')
    rna_sc_ad.write(snakemake.config['anndata']['sc_rna'])

    # CReate directory to mark completion
    import os
    os.makedirs(f"{str(snakemake.config['output'])}/gene_x_tf", exist_ok=True)
    with open(snakemake.output, 'a'):
        os.utime(snakemake.output, None)

if __name__ == "__main__":
    main()
