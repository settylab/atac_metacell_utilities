if __name__ == "__main__":
    import argparse
    desc = "Resize the widths of peaks"
    
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
    parser.add_argument(
        "--width",
        type=int,
        help="Width of final peak",
        default=150,
        metavar="width"
    )
    
    args = parser.parse_args()

import pyranges as pr
import pandas as pd
import scanpy as sc
import numpy as np
# If we make the `_pyranges_from_strings()` function in SEACells public, I can
#     I can use that instead of defining this function!
def pyranges_from_strings(pos_list):
    # Chromosome and positions
    chrom = pos_list.str.split(':').str.get(0)
    start = pd.Series(pos_list.str.split(':').str.get(1)).str.split('-').str.get(0)
    end = pd.Series(pos_list.str.split(':').str.get(1)).str.split('-').str.get(1)
    
    # Create ranges
    gr = pr.PyRanges(chromosomes=chrom, starts=start, ends=end)
    
    return gr


def make_peak_df(peaks_pr):
    
    peaks_df = pd.DataFrame()
    
    midpoint = peaks_pr.Start + ((peaks_pr.End - peaks_pr.Start) / 2).astype(int)
    start = midpoint - np.floor(args.width / 2).astype(int)
    end = midpoint + np.floor(args.width / 2).astype(int)
    
    # Positions
    peaks_df['chrom'] = peaks_pr.Chromosome
    peaks_df['chromStart'] = start
    peaks_df['chromEnd'] = end
    
    # summit
    peaks_df['summit']  = (args.width / 2).astype(int)

    # Score
    peaks_df['score'] = 1

    # Names
    peaks_df['name'] = peaks_df['chrom'].astype(str) + ':' + peaks_df['chromStart'].astype(str) + '-' + peaks_df['chromEnd'].astype(str)

    return peaks_df

def main(args):
    # Load data
    atac_ad = sc.read(args.atac)

    peaks_pr = pyranges_from_strings(atac_ad.var_names)
    
    peaks_df = make_peak_df(peaks_pr)
    peaks_df.to_csv(args.outdir + 'peaks.bed', sep='\t', index=None, header=True)
    
    
if __name__ == "__main__":
    main(args)
    
    