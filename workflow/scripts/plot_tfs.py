if __name__ == "__main__":
    import argparse
    desc = "Plots chromVAR scores for set of TFs"
    
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--sc_ad",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to scAnnData with embedding and annotations",
    )
    
    parser.add_argument(
        "--chromvar_outdir",
        metavar="directory",
        type=str,
        required=True,
        help="Path to chromVAR output directory",
    )
    
    parser.add_argument(
        "--tf_list",
        metavar="list",
        nargs="*",
        type=str,
        required=True,
        help="List of TFs to plot",
    )
    
    parser.add_argument(
        "--obs_annos",
        metavar="list",
        nargs="*",
        type=str,
        required=True,
        help="List of .obs keys to plot",
    )
    
    parser.add_argument(
        "--embedding",
        metavar="str",
        type=str,
        default='umap',
        help="embedding key in obsm to use",
    )
    
    parser.add_argument(
        "--vmax",
        metavar="int",
        type=int,
        default='3',
        help="vmax for plotting, vmin = -vmax",
    )
    
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Path to output directory",
        default="results/tf_plots/",
        metavar="directory"
    )

    args = parser.parse_args()

import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as plt


def main(args):
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    
    print('reading in chromVAR outputs...')
    cv_dev = pd.read_csv(f'{args.chromvar_outdir}/deviations.csv', index_col=0).T
    cv_dev_zs = pd.read_csv(f'{args.chromvar_outdir}/deviations_zs.csv', index_col=0).T 
    
    sc_ad = sc.read(args.sc_ad)
    
    print('creating AnnData...')
    cv_ad = sc.AnnData(cv_dev, obsm=sc_ad[cv_dev.index,:].obsm, obs=sc_ad[cv_dev.index,:].obs)
    cv_ad.layers['z-scored'] = cv_dev_zs.values
    
    if len(args.tf_list) != 0: 
        tfs_plot = [x for x in args.tf_list if x in cv_ad.var_names]
        not_valid = [x for x in args.tf_list if x not in cv_ad.var_names]
        
        if len(tfs_plot) != 0:
            print(f'Plotting the following TFs: {tfs_plot}')
            if len(not_valid) != 0:
                print(f'Following TFs not in chromVAR output: {not_valid}')            
            
            sc.pl.embedding(cv_ad, basis=args.embedding, color=tfs_plot, layer='z-scored', vmin=(-1 * args.vmax), vmax=args.vmax, 
                            cmap='RdBu_r', show=False, frameon=False)
            plt.savefig(f'{args.outdir}/tfs_zs.png', bbox_inches='tight', dpi=300)

            sc.pl.embedding(cv_ad, basis=args.embedding, color=tfs_plot, vmin=(-1 * args.vmax), vmax=args.vmax, 
                            cmap='RdBu_r', show=False, frameon=False)
            plt.savefig(f'{args.outdir}/tfs.png', bbox_inches='tight', dpi=300)        
        
        else:
            print('No valid TFs are in chromVAR output, supply new list')
    
    annos_plot = [x for x in args.obs_annos if x in sc_ad.obs.columns]
    if len(annos_plot) != 0:
        sc.pl.embedding(cv_ad, basis=args.embedding, color=annos_plot, frameon=False, show=False) 
        plt.savefig(f'{args.outdir}/annotated.png', bbox_inches='tight', dpi=300)
                                 

if __name__ == "__main__":
    main(args)