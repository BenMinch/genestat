import argparse
import itertools
import csv
from collections import defaultdict
from Bio import SeqIO
import os,sys,shutil
from glob import glob
import pandas as pd


def generate_tetranucleotides():
    """Generate all possible 4-mer combinations from A, C, T, G."""
    return ["".join(kmer) for kmer in itertools.product("ACTG", repeat=4)]

def count_tetranucleotides(genome_file):
    """Count tetranucleotide frequencies in a genome FASTA file."""
    tetranucleotide_counts = defaultdict(int)
    total_kmers = 0
    
    # Initialize dictionary with all possible 4-mers set to 0
    for kmer in generate_tetranucleotides():
        tetranucleotide_counts[kmer] = 0
    
    # Read genome sequence from FASTA
    with open(genome_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            sequence = str(record.seq).upper()
            for i in range(len(sequence) - 3):  # Sliding window of size 4
                kmer = sequence[i:i+4]
                if kmer in tetranucleotide_counts:
                    tetranucleotide_counts[kmer] += 1
                    total_kmers += 1
    
    # Normalize to proportions
    for kmer in tetranucleotide_counts:
        tetranucleotide_counts[kmer] = tetranucleotide_counts[kmer] / total_kmers if total_kmers > 0 else 0
    
    return tetranucleotide_counts

def write_output(output_file, tetranucleotide_frequencies):
    """Write tetranucleotide frequencies to a CSV file."""
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["4-mer", "Frequency"])
        for kmer, frequency in sorted(tetranucleotide_frequencies.items()):
            writer.writerow([kmer, frequency])

def merge_tetranucleotide_csvs(input_folder, output_file):
    """Merge multiple tetranucleotide frequency CSV files into one matrix."""
    
    # Get all CSV files in the folder
    csv_files = glob(os.path.join(input_folder, "*.csv"))
    
    if not csv_files:
        print("No CSV files found in the directory!")
        return

    # Initialize an empty dictionary to store dataframes
    dataframes = {}

    for file in csv_files:
        genome_name = os.path.basename(file).replace(".csv", "")  # Extract genome name
        df = pd.read_csv(file)

        # Ensure proper column names
        if df.columns.tolist() != ["4-mer", "Frequency"]:
            print(f"Skipping {file} - Incorrect column format")
            continue
        
        df.set_index("4-mer", inplace=True)  # Use tetranucleotide as index
        dataframes[genome_name] = df["Frequency"]

    # Merge all dataframes on the "4-mer" index
    merged_df = pd.DataFrame(dataframes).fillna(0)  # Fill missing values with 0

    # Reset index to move "4-mer" back to a column
    merged_df.reset_index(inplace=True)

    # Save to CSV
    #transpose the dataframe to have genomes as rows and tetranucleotides as columns
    merged_df = merged_df.T
    merged_df.columns = merged_df.iloc[0]  # Set the first row as header
    merged_df = merged_df[1:]  # Remove the first row
    merged_df.index.name = "Genome"  # Set index name
    merged_df.reset_index(inplace=True)  # Reset index to make "Genome" a column
    merged_df.columns.name = None  # Remove the index name
    merged_df = merged_df.rename_axis(None, axis=1)  # Remove the index name
    merged_df.to_csv(output_file, index=False)
    print(f"Merged file saved as: {output_file}")



def main():
    parser = argparse.ArgumentParser(description="Calculate tetranucleotide frequencies from a genome FASTA file.")
    parser.add_argument('-i',"--genome", help="Input genome file in FASTA format or directory of genomes")
    parser.add_argument('-o',"--output", help="Output CSV file to store 4-mer frequencies")
    parser.add_argument('-b',"--batch",help="Batch mode",action="store_true",default=False,required=False)
    args = parser.parse_args()
    
    if args.batch:
        if not os.path.exists(args.output):
            os.makedirs(args.output)
        print("Batch mode is enabled")
        for genome in os.listdir(args.genome):
            path=os.path.join(args.genome,genome)
            tetranucleotide_frequencies = count_tetranucleotides(path)
            output_file = os.path.join(args.output, genome + "TNF.csv")
            write_output(output_file, tetranucleotide_frequencies)
        #Combine all the files into one
        merge_tetranucleotide_csvs(args.output, "all_genomes_tetranucleotides_"+args.output+".csv")
    else:
        print("Processing the genome file: ",args.genome)
        print("Output file: ",args.output)
        tetranucleotide_frequencies = count_tetranucleotides(args.genome)
        write_output(args.output, tetranucleotide_frequencies)

    



if __name__ == "__main__":
    main()
