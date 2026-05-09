"""Split a multi-record FASTA into N roughly-equal chunks for parallel FIMO.

Round-robin assignment by record index (not by base count) keeps chunk sizes
balanced for the typical case where all sequences are the same span (FIMO
input from seq_gl.py / all_seqs.fa is uniform width). Two-pass implementation
keeps memory bounded.
"""

import argparse
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input_fa", help="Input FASTA")
    p.add_argument(
        "out_prefix",
        help="Output prefix; chunks written as {prefix}{i}.fa for i in 0..n-1",
    )
    p.add_argument("n", type=int, help="Number of chunks")
    return p.parse_args()


def split_fasta(input_fa, out_prefix, n):
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    out_dir = os.path.dirname(out_prefix) or "."
    os.makedirs(out_dir, exist_ok=True)

    handles = [open(f"{out_prefix}{i}.fa", "w") for i in range(n)]
    counts = [0] * n
    try:
        record_idx = -1
        current = None
        with open(input_fa, "r") as f:
            for line in f:
                if line.startswith(">"):
                    record_idx += 1
                    current = record_idx % n
                    handles[current].write(line)
                    counts[current] += 1
                else:
                    if current is None:
                        raise ValueError(
                            f"FASTA {input_fa} starts with non-header line"
                        )
                    handles[current].write(line)
    finally:
        for h in handles:
            h.close()

    total = sum(counts)
    print(f"split {total} records into {n} chunks: {counts}", file=sys.stderr)


if __name__ == "__main__":
    args = parse_args()
    split_fasta(args.input_fa, args.out_prefix, args.n)
