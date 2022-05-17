#! /usr/bin/env python
import os
import argparse

desc = """
Make a MEME-format motif file from a CIS-BP data dump.

Example: make_meme.py \
    --motif-table /fh/fast/setty_m/user/dotto/prime/data/cisbp_data/TF_Information.txt \
    --pwms-dir /fh/fast/setty_m/user/dotto/prime/data/cisbp_data/pwms_all_motifs \
    --base-url "http://cisbp.ccbr.utoronto.ca/TFreport.php?versionNumber=2.00&searchTF="
    data/cis-bp-tf-information.meme
"""
parser = argparse.ArgumentParser(description=desc)
parser.add_argument(
    "out",
    metavar="file",
    type=str,
    nargs="?",
    help="outputfile",
)

parser.add_argument(
    "--motif-table",
    help="CIS-BP TF_Information table",
    type=str,
    default=None,
    metavar="file",
)
parser.add_argument(
    "--sep",
    help="Sperator used in --motif-table",
    type=str,
    default="\t",
    metavar="URL",
)
parser.add_argument(
    "--pwms-dir",
    help="Directory with CIS-BP pwms files.",
    type=str,
    default=None,
    metavar="directory",
)
parser.add_argument(
    "--base-url",
    help="URL to use in the output file",
    type=str,
    default=None,
    metavar="URL",
)

args = parser.parse_args()


def read_pwms(file):
    first = True
    lines = list()
    with open(file) as fl:
        for line in fl.readlines():
            if first:
                first = False
                continue
            data = [
                float(p) for i, p in enumerate(line.strip().split(args.sep)) if i > 0
            ]
            total = sum(data)
            data = [str(p/total) for p in data]
            lines.append(" ".join(data))
    return lines


meme_start = """MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 T 0.25000 

"""


if __name__ == "__main__":
    header = None
    motif_count = 0
    print(f"Writing {args.out}")
    with open(args.motif_table, "r") as tfl:
        with open(args.out, "w") as outfl:
            outfl.write(meme_start)
            for raw_line in tfl.readlines():
                line = raw_line.strip()
                if len(line) == 0:
                    continue
                if header is None:
                    header = line.split(args.sep)
                    for idx, col in enumerate(header):
                        if col == "Motif_ID":
                            motif_id_idx = idx
                        elif col == "TF_Name":
                            tf_name_idx = idx
                        elif col == "TF_ID":
                            tf_id_idx = idx
                    idx_tupel = (motif_id_idx, tf_name_idx, tf_id_idx)
                    continue
                new_lines = list()
                data = line.split(args.sep)
                motif_id, tf_name, tf_id = [data[i] for i in idx_tupel]
                if motif_id == ".":
                    continue
                new_lines.append(f"MOTIF {motif_id} {tf_name}")
                new_lines.append("")
                pwms_file = os.path.join(args.pwms_dir, motif_id + ".txt")
                probabilities = read_pwms(pwms_file)
                width = len(probabilities)
                if width == 0:
                    continue
                new_lines.append(
                    f"letter-probability matrix: alength= 4 w= {width} "
                    "nsites= 1 E= 0"
                )
                new_lines += probabilities + [""]
                if args.base_url is not None:
                    new_lines.append(f"URL {args.base_url}{tf_id}")
                    new_lines.append("")
                new_lines.append("")
                outfl.write(os.linesep.join(new_lines))
                motif_count += 1
        print(f'{motif_count} motifs written.')
