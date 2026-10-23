#if __name__ == "__main__":
#    import argparse
#    desc = "Compute correlations between RNA gene expression and ATAC accessibility"
#
#    parser = argparse.ArgumentParser(
#        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
#    )
#
#    parser.add_argument(
#        "--atac",
#        metavar="AnnData",
#        type=str,
#        required=True,
#        help="Path to ATAC MC AnnData, common obs with RNA"
#    )
#
#    parser.add_argument(
#        "--rna",
#        metavar="AnnData",
#        type=str,
#        required=True,
#        help="Path to RNA MC AnnData, common obs with RNA"
#    )
#
#    parser.add_argument(
#        "-o",
#        "--outdir",
#        type=str,
#        help="Path to output directory",
#        default="results/gp_corrs/",
#        metavar="directory"
#    )
#
#    parser.add_argument(
#        "--gtf_file",
#        type=str,
#        help="Number of jobs for computing GP correlations",
#        default="data/hg38.gtf",
#        metavar="directory"
#    )
#
#    parser.add_argument(
#        "--n_jobs",
#        type=int,
#        help="Number of jobs for computing GP correlations",
#        default=1,
#        metavar="int"
#    )
#
#    parser.add_argument(
#        "--transcript_span",
#        type=int,
#        help="num of base pairs to add on longest transcript",
#        default=100000,
#        metavar="int"
#    )
#
#    parser.add_argument(
#        "--test_set",
#        help="Test on a subset of genes, n=--n_genes",
#        action="store_true"
#    )
#
#    parser.add_argument(
#        "--n_genes",
#        type=int,
#        help="num of base pairs to add on longest transcript",
#        default=20,
#        metavar="int"
#    )
#
#    args = parser.parse_args()

import scanpy as sc
import numpy as np
import SEACells
import pandas as pd
from tqdm.auto import tqdm
import os
import warnings

def main():
    
    out_path = str(snakemake.output[0])
    out_dir = os.path.dirname(out_path)
    print(f"Saving output to {out_dir}")
    np.seterr(divide = 'ignore')
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    print('Loading data...')
    atac_ad = sc.read(snakemake.config["anndata"]["atac"])
    rna_ad = sc.read(snakemake.config["anndata"]['rna'])
    if snakemake.params['test_set']:
        gene_set = rna_ad.var_names[:snakemake.params['n_genes']]
    else:
        gene_set = rna_ad.var_names

    # Compute Gene-Peak Correlation scores
    print('Computing gene peak correlations...')
    gp_corr = SEACells.genescores.get_gene_peak_correlations(atac_ad, rna_ad,
                                                             path_to_gtf=str(snakemake.input['gtf_file']),
                                                             span=snakemake.params['transcript_span'], n_jobs=snakemake.params['n_jobs'],
                                                             gene_set=gene_set)

    import pickle
    with open(out_dir + '/gp_corr.pickle', 'wb') as handle:
        pickle.dump(gp_corr, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Save Gene-peak correlations to Anndata
    gp_corr_df = pd.DataFrame()
    for gene in tqdm(gp_corr.keys()):
        if type(gp_corr[gene]) != int:
            gp_corr[gene]["gene"] = pd.Series(gene, gp_corr[gene].index)
            gp_corr_df = pd.concat([gp_corr_df, gp_corr[gene]])
    atac_ad.uns["gp_corrs"] = gp_corr_df

    # Write anndata file
    print('Saving anndata file ...')
    atac_ad.write(snakemake.config["anndata"]["atac"])

    
    os.makedirs(out_dir, exist_ok=True)
    with open(snakemake.output[0], 'a'):
        os.utime(snakemake.output[0], None)


if __name__ == "__main__":
    main()
