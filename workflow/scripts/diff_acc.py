import scanpy as sc
import pandas as pd
import os


import subprocess


def main():
    # Arguments
    to_compare = snakemake.params["to_compare"].split('/')
    data_dir = snakemake.params["out_dir"][0]

    # Load anndata
    atac_ad = sc.read(snakemake.config["anndata"]["atac"])
    #command = 'Rscript'
    #path_to_script = str(snakemake.input['script'])

    # Run the script
    #subprocess.check_call([command, path_to_script, str(data_dir), str(snakemake.params['to_compare']), str(snakemake.params['cell_type_obs'])], shell=False)

    # Load diff results and save to anndata
    for i in range(len(to_compare)):
        # Cell types
        to_compare[i] = to_compare[i].split(',')

        # Diff results
        df = pd.read_csv(data_dir + f"/{to_compare[i][0]}_{to_compare[i][1]}_diff_acc.tsv", sep="\t")
        df.index = df['feature']

        # Add to anndata
        atac_ad.varm[f'{to_compare[i][0]}_{to_compare[i][1]}_diff_acc'] = df.loc[atac_ad.var_names, :]

    # Save anndata
    atac_ad.write(snakemake.config["anndata"]["atac"])

    
        
    out_path = str(snakemake.output[0])
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'a'):
        os.utime(out_path, None)

if __name__ == "__main__":
    main()
