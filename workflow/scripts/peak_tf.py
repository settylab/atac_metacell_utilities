if __name__ == "__main__":
    import argparse
    desc = "Create a peaks X TF AnnData with FIMO scores"
    
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "atac",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to ATAC MC AnnData, common obs with RNA"
    )
    
    parser.add_argument(
        "--peak_file",
        metavar="file",
        type=str,
        required=True,
        help="Path to peaks.bed file",
    )
    
    parser.add_argument(
        "--fimo_res",
        metavar="directory",
        type=str,
        default="results/fimo_result/",
        help="Path to FIMO output directory",
    )
    
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Path to output directory",
        default="results/",
        metavar="directory"
    )
    args = parser.parse_args()

from tqdm.auto import tqdm  
from scipy.sparse import csr_matrix
import scanpy as sc
import pandas as pd
import subprocess

def build_peak_tf(fimo_scores, peaks_df, atac_ad):
    # Motif information 
    motifs = pd.Series()
    motif_index = 0
    rec_index = 0

    # Peak index
    peak_index = pd.Series(range(len(peaks_df.name)), index=peaks_df.name)

    # Initialize Values
    num_records = int(subprocess.run(['wc', '-l', fimo_res_path], stdout=subprocess.PIPE).stdout.decode().split(' ')[0]) - 5
    
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
    pxtf_scores = csr_matrix((values, (x, y)), (atac_ad.shape[1], motif_index))
    pxtf_ad = sc.AnnData(pxtf_scores)
    pxtf_ad.var_names = motifs.index
    
    return pxtf_ad


def main(args):
    # Load data
    atac_ad = sc.read(args.atac)

    peaks_df = pd.read_csv(args.peak_file, sep='\t')
    fimo_scores = args.fimo_res + 'fimo.tsv'
    
    # Build peak x TF AnnData
    pxtf_ad = build_peak_tf(fimo_scores, peaks_df, atac_ad)

    pxtf_ad.obs_names = peaks_df.name
    
    # annotated the full size peak
    pxtf_ad.obs['unsized_peak'] = atac_ad.var_names
    
    pxtf_ad.write(args.outdir + 'peak_x_tf.h5ad')

    
if __name__ == "__main__":
    main(args)