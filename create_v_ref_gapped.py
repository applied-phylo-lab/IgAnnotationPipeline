import argparse
import pandas as pd
import os

def main(input_file, output_dir, locus):
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t', dtype=str)

    # Check required columns exist
    if 'v_sequence_alignment' not in df.columns or 'sequence_id' not in df.columns:
        raise ValueError("Input file must contain 'v_sequence_alignment' and 'sequence_id' columns.")

    # Drop rows where v_sequence_alignment is missing
    df = df.dropna(subset=['v_sequence_alignment', 'sequence_id'])

    # Prepare output path
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{locus}V_gapped.fasta")

    # Write to FASTA
    with open(output_file, "w") as fasta_out:
        for _, row in df.iterrows():
            fasta_out.write(f">{row['sequence_id']}\n{row['v_sequence_alignment'].upper()}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract V-gapped sequences to FASTA.")
    parser.add_argument("-i", "--input", required=True, help="Input .tsv file from VQUEST")
    parser.add_argument("-o", "--output", required=True, help="Output directory for FASTA")
    parser.add_argument("-l", "--locus", required=True, help="Locus name (e.g. IGH, IGK)")

    args = parser.parse_args()

    main(args.input, args.output, args.locus)
