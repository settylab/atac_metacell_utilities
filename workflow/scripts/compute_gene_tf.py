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
        "--sc_atac",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to ATAC single-cell AnnData, common obs with RNA",
    )

    parser.add_argument(
        "--sc_rna",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to RNA single-cell AnnData, common obs with ATAC",
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
import scanpy as sc
from scipy.sparse import csr_matrix


def main(args):
    # Load data
    print('Loading data')
    atac_sc_ad = sc.read(args.sc_atac)
    atac_mc_ad = sc.read(args.atac)
    rna_sc_ad = sc.read(args.sc_rna)

    # Assemble necessary matrices
    print('Computing gene TF relationships')

    # Gene X Peak matrix
    gp = atac_mc_ad.uns['gp_corrs'].copy()
    gp = gp.loc[(gp['cor'] > args.min_corr) & (gp['pval'] < args.max_pval), :]
    peak_counts = gp.groupby('gene').apply(len)
    use_genes = peak_counts[peak_counts >= args.min_peaks].index
    gp.loc[gp['gene'].isin(use_genes), :]

    # Convert to categoricals
    gp['peaks'] = gp.index.values
    gp['peaks'] = pd.Categorical(gp['peaks'], atac_sc_ad.var_names)
    gp['gene'] = pd.Categorical(gp['gene'], rna_sc_ad.var_names)
    
    # Binarize gene peak correlations
    gp['cor'] = (gp['cor'] > 0).astype(int)
    
    # Construct sparse matrix
    data = gp['cor'].values
    row_ind = gp['gene'].values.codes
    col_ind = gp['peaks'].values.codes
    M = len(gp['gene'].values.categories)
    N = len(gp['peaks'].values.categories)
    gene_X_peak = csr_matrix((data, (row_ind, col_ind)), [M, N])

    # Peak X TF matrix
    peak_X_tf = atac_sc_ad.varm['InSilicoChip']

    # Gene X TF matrix
    gene_X_tf = gene_X_peak.dot(peak_X_tf)
    
    rna_sc_ad.varm['geneXTF'] = pd.DataFrame(gene_X_tf.todense(),
        index=rna_sc_ad.var_names, columns=atac_sc_ad.uns['InSilicoChipColumns'])
    
    # Saving resuls
    print('Saving results')
    rna_sc_ad.write(args.sc_rna)

    # CReate directory to mark completion
    import os
    os.makedirs(str(args.outdir), exist_ok=True)

if __name__ == "__main__":
    main(args)
