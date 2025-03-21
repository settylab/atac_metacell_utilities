import os
import sys
import argparse

def split_meme_file(input_file, output_dir):
    """
    Splits a combined MEME file into individual .meme files.

    Parameters:
        input_file (str): Path to the combined MEME file.
        output_dir (str): Directory to save the individual .meme files.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(input_file, 'r') as f:
        lines = f.readlines()

    header = []
    motifs = []
    current_motif = []
    in_motif_section = False

    for line in lines:
        # Store the header until the start of the first motif
        if line.startswith('MOTIF'):
            in_motif_section = True
            if current_motif:
                motifs.append(current_motif)
                current_motif = []
        if in_motif_section:
            current_motif.append(line)
        else:
            header.append(line)

    # Add the last motif
    if current_motif:
        motifs.append(current_motif)

    # Write each motif to its own file
    for i, motif_lines in enumerate(motifs):
        motif_name = motif_lines[0].split()[1]
        output_file = os.path.join(output_dir, f"{motif_name}.meme")
        with open(output_file, 'w') as f:
            f.writelines(header)
            f.writelines(motif_lines)

    print(f"Successfully split {len(motifs)} motifs into individual .meme files in {output_dir}")


if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Split a combined MEME file into individual .meme files.")
    parser.add_argument("input_file", help="Path to the combined MEME file.")
    parser.add_argument("output_dir", help="Directory to save individual .meme files.")
    
    # Parse arguments
    args = parser.parse_args()

    # Run the split function
    split_meme_file(args.input_file, args.output_dir)

