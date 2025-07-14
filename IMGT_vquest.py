import argparse
import subprocess
import re
import sys

# Dictionary to map species codes to IMGT species names
SPECIES_MAP = {
    "human": "human",
    "mLemCat1_alternate_igdetective": "lemur",
    "mLemCat1_primary_igdetective": "lemur",
    "mPonAbe1_hap1": "pongo",
    "mPonPyg2_hap1": "pongopygmaeus",
    "mGorGor1_mat": "gorilla",
    "mNycCou1_primary_igdetective": "lemur",
    "mNycCou1_alternate_igdetective": "lemur",
    "mPanPan1_mat": "human",
    "Macaca_mulatta": "rhesus-monkey"
}

def extract_species(species_param):
    """
    Extracts the species identifier from the full parameter by removing locus and gene type.
    
    Examples:
    'mLemCat1_alternate_igdetective_IGH_V' -> 'mLemCat1_alternate_igdetective'
    'mLemCat1_alternate_igdetective_V' -> 'mLemCat1_alternate_igdetective'
    """
    return re.sub(r'(_[A-Z]+)?_[A-Z]+$', '', species_param)


def main():
    parser = argparse.ArgumentParser(description="Run vquest with the correct species mapping.")
    parser.add_argument("-f", "--fasta", required=True, help="Path to input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-l", "--locus", required=True, help="Locus type")
    parser.add_argument("-s", "--species", required=True, help="Species identifier")

    args = parser.parse_args()

    # Extract base species name
    species_base = extract_species(args.species)

    # Get the IMGT species name
    species_matched = SPECIES_MAP.get(species_base)
    if species_matched is None:
        print(f"Error: No IMGT species mapping found for '{species_base}'")
        exit(1)

    # Construct and run the vquest command
    command = [
        "vquest",
        "--species", species_matched,
        "--receptorOrLocusType", args.locus,
        "--fileSequences", args.fasta,
        "-o", args.output
    ]

    print(f"Running command: {' '.join(command)}")
    try: 
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
    # Log error, create placeholder output
        file= args.output+"vquest_airr.tsv"
        with open(file, "w") as f:
            f.write("")
        sys.exit(0)

if __name__ == "__main__":
    main()
