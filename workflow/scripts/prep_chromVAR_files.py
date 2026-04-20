import pandas as pd
import scipy.io
from scipy.sparse import csr_matrix
import scanpy as sc
import os
from snakemake.script import snakemake

def filt_and_binarize(ins_chip_mat, min_chip_score, min_peak_hits):
    print(f'Filtering for TFs with insilico chip score > {min_chip_score} in at least {min_peak_hits} peaks')
    ins_chip_mat[ins_chip_mat <= min_chip_score] = 0.0
    ins_chip_bin = (ins_chip_mat != 0.0).astype(int)
    filter_bin = ins_chip_bin.loc[:, ins_chip_bin.sum() >= min_peak_hits]

    return filter_bin


def main():
    # Load ATAC data and reconstruct in-silico chip matrix
    atac_sc_ad = sc.read(snakemake.config["anndata"]["sc_atac"])

    # insilio Chip results
    insc_df = pd.DataFrame(atac_sc_ad.varm['InSilicoChip'].todense(),
                           columns=atac_sc_ad.uns['InSilicoChipColumns'], index=atac_sc_ad.var_names)

    binary_mat = filt_and_binarize(insc_df, snakemake.params['min_chip'], snakemake.params['min_peak_hits'])
    print(f'Filtered {insc_df.shape[1]} down to {binary_mat.shape[1]} TFs for chromVAR')

    # Binary for chromvar
    print(f'writing binary in silico ChIP matrix and TF names...')
    os.makedirs(f"{snakemake.config["output"]}chromvar_input/", exist_ok=True)
    # export binary matrix for chromVAR
    scipy.io.mmwrite(snakemake.output['bin_ins_chip'],
                     csr_matrix(binary_mat.values))
    # export tf names
    pd.Series(binary_mat.columns).to_csv(snakemake.output['chromvar_tf_names']
                                         , index=False, header=['tf_name'])

    # Peak counts
    scipy.io.mmwrite(snakemake.output['sc_atac_counts'], atac_sc_ad.X)
    atac_sc_ad.obs['nFrags'].to_csv(
        snakemake.output['sc_nfrags'], index=True, header=['nFrags'])


if __name__ == "__main__":
    main()
