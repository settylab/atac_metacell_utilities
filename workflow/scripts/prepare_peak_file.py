import pandas as pd
import scanpy as sc
import mudata as md


def make_peak_df(peaks_list):

    # Information about peaks
    chrom = peaks_list.str.split(':').str.get(0)
    start = pd.Series(peaks_list.str.split(
        ':').str.get(1)).str.split('-').str.get(0)
    end = pd.Series(peaks_list.str.split(
        ':').str.get(1)).str.split('-').str.get(1)

    # Positions
    peaks_df = pd.DataFrame()
    peaks_df['chrom'] = chrom
    peaks_df['chromStart'] = start.astype(int)
    peaks_df['chromEnd'] = end.astype(int)

    # Summit
    peaks_df['summit'] = (
        (peaks_df['chromEnd'] - peaks_df['chromStart']) / 2).astype(int)

    # Score
    peaks_df['score'] = 1

    # Names
    peaks_df['name'] = peaks_list.values

    return peaks_df


def main():
    # Load data
    atac_ad = md.read(snakemake.input.atac)
    peaks_df = make_peak_df(atac_ad.var_names)
    peaks_df.to_csv(snakemake.params.out_dir + 'peaks.bed',
                    sep='\t', index=None, header=False)


if __name__ == "__main__":
    main()
