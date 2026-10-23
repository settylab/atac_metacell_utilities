from tqdm.auto import tqdm
from scipy.sparse import csr_matrix
import scanpy as sc
import pandas as pd
import numpy as np
import subprocess

def build_peak_tf(fimo_scores, peaks_df):
    # Motif information
    motifs = pd.Series(dtype=object)
    motif_index = 0
    rec_index = 0
    
    # Peak index
    peak_index = pd.Series(range(len(peaks_df.name)), index=peaks_df.name)

    # Initialize Values
    num_records = int(subprocess.run(
        ['wc', '-l', fimo_scores], stdout=subprocess.PIPE).stdout.decode().split(' ')[0]) - 5
    
    x = np.zeros(num_records)
    y = np.zeros(num_records)
    values = np.zeros(num_records)

    # Read file

    with open(fimo_scores, 'r') as f:
        for line in tqdm(f):
            # Skip first line
            split = line.split('\t')
            if split[0] == 'motif_id':
                continue

            if len(split) == 1:
                break

            # Update motifs if necessary
            if split[1] not in motifs:
                motifs[split[1]] = motif_index
                motif_index += 1

            # Update record
            x[rec_index] = peak_index[split[2]]
            y[rec_index] = motifs[split[1]]
            values[rec_index] = float(split[6])
            rec_index += 1

    # Create AnnData
    pxtf_scores = csr_matrix(
        (values, (x, y)), (len(peaks_df.name), motif_index))
    pxtf_ad = sc.AnnData(pxtf_scores)
    pxtf_ad.var_names = motifs.index
    pxtf_ad.obs_names = peaks_df['name']

    return pxtf_ad


def main():
    print(snakemake.output)
    peaks_file = str(snakemake.input["peaks"])
    peaks_df = pd.read_csv(peaks_file, sep='\t', header=None)
    peaks_df.columns = ['chrom', 'chromStart',
                        'chromEnd', 'summit', 'score', 'name']
    fimo_scores = str(snakemake.input['fimo_file'])

    # Build peak x TF AnnData
    print('Building peak X TF matrix')
    pxtf_ad = build_peak_tf(fimo_scores, peaks_df)

    # Load and save in anndata
    print('Saving results to anndata')
    atac_sc_ad = sc.read(snakemake.config['anndata']["sc_atac"])
    atac_sc_ad.varm['FIMO'] = pxtf_ad[atac_sc_ad.var_names, :].X
    atac_sc_ad.uns['FIMOColumns'] = np.array(pxtf_ad.var_names)
    atac_sc_ad.write(snakemake.config['anndata']["sc_atac"])

    # CReate directory to mark completion
    import os
    os.makedirs(f"{str(snakemake.config['output'])}peak_tf/", exist_ok=True)
    out_file = str(snakemake.output)
    with open(out_file, 'a'):
        os.utime(out_file, None)


if __name__ == "__main__":
    main()
