#if __name__ == "__main__":
#    import argparse
#    desc = "Adds Chromvar results to ATAC annData."
#
#    parser = argparse.ArgumentParser(
#        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
#    )
#
#    parser.add_argument(
#        "--sc_rna",
#        metavar="AnnData",
#        type=str,
#        required=True,
#        help="Path to RNA single-cell AnnData",
#    )
#    parser.add_argument(
#        "--input",
#        type=str,
#        required=True,
#        help="Path to ChromVAR results for input.",
#    )
#
#    args = parser.parse_args()
#
import scanpy as sc
import pandas as pd
#import subprocess

def main():
    
    #print("Running run_chromvar.R")
    #command = 'Rscript'
    #path_to_script = snakemake.input['script']

    # Run the script
    # {input.script} {input.peak_file} {params.chromvar_indir} {output.chromvar_outdir} {params.genome}
    #subprocess.check_call([command, path_to_script, snakemake.input['peak_file'], snakemake.params['chromvar_indir'], snakemake.params['chromvar_outdir'], snakemake.params['genome'], shell=False)
    print('Loading and then saving chromVAR results to anndata..')
    rna_ad = sc.read_h5ad(snakemake.input['sc_rna'])
    deviations = pd.read_csv(args.input + "/deviations.csv", index_col=0).T
    deviations.index = rna_ad.obs_names
    rna_ad.obsm["chromVAR_deviations"] = deviations
    rna_ad.write(args.sc_rna)


if __name__ == "__main__":
    main()
