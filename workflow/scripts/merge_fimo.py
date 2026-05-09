"""Concatenate per-chunk FIMO `fimo.tsv` files and recompute global q-values.

FIMO's TSV layout: motif_id, motif_alt_id, sequence_name, start, stop,
strand, score, p-value, q-value, matched_sequence -- with one header
line at the top and a few `#`-prefixed comment lines at the bottom
(version banner, command line). Under scatter, per-chunk q-values are
incorrect because q-value is a global FDR statistic. We run FIMO with
--no-qvalue per chunk (fast; col 9 emitted empty), then recompute
global q-values here via Benjamini-Hochberg over the concatenated
p-value column. Output schema is byte-identical to a single-job
unflagged FIMO TSV; downstream peak_tf.py is unaffected.
"""

import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output_tsv", help="Merged output path")
    p.add_argument("chunk_tsvs", nargs="+", help="Per-chunk fimo.tsv paths")
    return p.parse_args()


def bh_qvalues(pvals):
    """Benjamini-Hochberg FDR. Returns q-values in original input order."""
    n = len(pvals)
    if n == 0:
        return []
    # Sort by p ascending; track original indices.
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    # First pass: q_i = p_i * n / rank_i (rank starts at 1).
    for rank_minus_1, idx in enumerate(order):
        q[idx] = pvals[idx] * n / (rank_minus_1 + 1)
    # Enforce monotonicity from the largest p downward, and clip to 1.
    prev = 1.0
    for idx in reversed(order):
        prev = min(prev, q[idx])
        q[idx] = prev
    return q


def merge(output_tsv, chunk_tsvs):
    header = None
    rows = []  # list of (cols list,) so we can rewrite col 8 (q-value)
    for chunk_idx, path in enumerate(chunk_tsvs):
        with open(path, "r") as f:
            for line in f:
                stripped = line.rstrip("\n")
                if not stripped:
                    continue
                if stripped.startswith("#"):
                    # Trailing version/command comments -- drop entirely.
                    continue
                if stripped.startswith("motif_id\t"):
                    if header is None:
                        header = stripped
                    continue
                rows.append(stripped.split("\t"))

    if header is None:
        # All chunks empty -- emit just the standard header for downstream
        # parsers that key off it.
        header = "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence"

    # Recompute global q-values from col index 7 (p-value).
    if rows:
        try:
            pvals = [float(r[7]) for r in rows]
        except (IndexError, ValueError) as e:
            raise RuntimeError(
                f"Could not parse p-value (column 8) from FIMO chunks: {e}"
            )
        qvals = bh_qvalues(pvals)
        for r, q in zip(rows, qvals):
            # Pad short rows defensively (shouldn't happen with default FIMO).
            while len(r) < 10:
                r.append("")
            r[8] = f"{q:.3g}"

    with open(output_tsv, "w") as out:
        out.write(header + "\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    print(
        f"merged {len(chunk_tsvs)} chunks -> {output_tsv} "
        f"({len(rows)} hits; q-values recomputed via global BH-FDR)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    args = parse_args()
    merge(args.output_tsv, args.chunk_tsvs)
