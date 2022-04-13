if __name__ == "__main__":
    import argparse
    desc = "Prepares peak file for downstream MOTIF analysis"
    
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--atac",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to ATAC AnnData.",
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

import pandas as pd
import scanpy as sc


def make_peak_df(peaks_list):

    # Information about peaks
    chrom = peaks_list.str.split(':').str.get(0)
    start = pd.Series(peaks_list.str.split(':').str.get(1)).str.split('-').str.get(0)
    end = pd.Series(peaks_list.str.split(':').str.get(1)).str.split('-').str.get(1)

    # Positions
    peaks_df = pd.DataFrame()
    peaks_df['chrom'] = chrom
    peaks_df['chromStart'] = start.astype(int)
    peaks_df['chromEnd'] = end.astype(int)

    # Summit
    peaks_df['summit'] = ((peaks_df['chromEnd'] - peaks_df['chromStart']) / 2).astype(int)

    # Score
    peaks_df['score'] = 1

    # Names
    peaks_df['name'] = peaks_list.values

    return peaks_df

def main(args):
    # Load data
    atac_ad = sc.read(args.atac)

    peaks_df = make_peak_df(atac_ad.var_names)
    peaks_df.to_csv(args.outdir + 'peaks.bed', sep='\t', index=None, header=True)
    
    
if __name__ == "__main__":
    main(args)
    
    
