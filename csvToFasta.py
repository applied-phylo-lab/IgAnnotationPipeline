import csv
import os
import sys
from collections import defaultdict
import pandas as pd

def main(input_csv, output_path,mapping,locus):
    # Read the Latin to Common name mapping file
    latin_to_common_file = mapping
    latin_to_common_mapping = {}

    with open(latin_to_common_file, mode="r") as mapping_file:
        mapping_reader = csv.DictReader(mapping_file)
        for row in mapping_reader:
            latin_to_common_mapping[row['FilePrefix']] = {
                "LatinName": row['LatinName'],
                "source": row['source']
            }

    # Determine the FilePrefix used in the input file name
    file_basename = os.path.basename(output_path)
    matched_prefix = None

    for prefix in latin_to_common_mapping:
        if prefix in file_basename:
            matched_prefix = prefix
            break


    if not matched_prefix:
        raise ValueError(f"No matching FilePrefix found in {latin_to_common_file} for file '{file_basename}'")

    latin_name = latin_to_common_mapping[matched_prefix]['LatinName']
    source = latin_to_common_mapping[matched_prefix]['source']

   # Detect delimiter and read the input file
    with open(input_csv, mode="r") as file:
        first_line = file.readline()
        delimiter = "\t" if "\t" in first_line else ","

    df = pd.read_csv(input_csv, delimiter=delimiter)
    df = df[df.columns].loc[df[df.columns[0]] != df.columns[0]]  # Remove duplicate header rows


    # Ensure required columns exist
    required_cols = ['Sequence', 'GeneType', 'Contig', 'Pos', 'Productive']
    for col in required_cols:
        if col not in df.columns:
            raise KeyError(f"The input CSV does not contain a required column: '{col}'")

    # Remove duplicate sequences
    df = df.drop_duplicates(subset=['Sequence'])

    v_df = df[df['GeneType'] == 'V']
    if v_df.empty:
        raise ValueError("No rows with GeneType == 'V' found in the input CSV.")
    contig_counts = v_df['Contig'].value_counts()
    main_contig = contig_counts.idxmax()

    # Filter entire DataFrame to only include rows from this contig
    df = df[df['Contig'] == main_contig]

    # Organize data by GeneType
    gene_type_sequences = defaultdict(list)

    for _, row in df.iterrows():
        row['source'] = source
        row['LatinName'] = latin_name

        gene_type = row['GeneType']
        header = f"{source}.{latin_name}.{row['Pos']}.{row['Contig']}.{gene_type}.{row['Productive']}"
        sequence = row['Sequence'].upper()

        gene_type_sequences[gene_type].append((header, sequence))

    # Write the sequences to separate FASTA files
    os.makedirs(output_path, exist_ok=True)

    for gene_type, sequences in gene_type_sequences.items():
        fasta_file_path = os.path.join(output_path, f"{locus}{gene_type}.fasta")
        with open(fasta_file_path, mode="w") as fasta_file:
            for header, sequence in sequences:
                fasta_file.write(f">{header}\n{sequence}\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python csvToFasta.py <input_csv> <output_path> <mapping_csv> <locus>")
        sys.exit(1)

    input_csv = sys.argv[1]
    output_path = sys.argv[2]
    mapping = sys.argv[3]
    locus = sys.argv[4]

    main(input_csv, output_path,mapping,locus)
