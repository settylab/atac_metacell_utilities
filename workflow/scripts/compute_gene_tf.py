if __name__ == "__main__":
    import argparse
    desc = "Creates three versions of gene x TF matrices"
    
    parser = argparse.ArgumentParser(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "atac",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to ATAC MC AnnData.",
    )
    
    parser.add_argument(
        "rna",
        metavar="AnnData",
        type=str,
        required=True,
        help="Path to RNA MC AnnData.",
    )
    
    parser.add_argument(
        "--peak_tf",
        metavar="file",
        type=str,
        required=True,
        help="Path to peak x TF AnnData",
    )
    
    parser.add_argument(
        "--gp_corr",
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
