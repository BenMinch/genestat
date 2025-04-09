import os
import re
import csv
import argparse

# Function to find internal stop codons in a gene sequence
def count_internal_stop_codons(sequence):
    stop_codons = ["TAA", "TAG", "TGA"]
    internal_stop_codons = 0
    sequence_length = len(sequence)

    # Check for stop codons within the sequence, excluding first and last codons
    for i in range(1, sequence_length - 1, 3):  # Step by 3 to check codons
        codon = sequence[i:i+3]
        if codon in stop_codons:
            internal_stop_codons += 1

    return internal_stop_codons

# Function to process each gene file in the folder
def process_genome_folder(folder_path):
    genome_stop_codons = []

    # Loop through all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".fasta") or filename.endswith(".fa")or filename.endswith(".fna"):  # Only process FASTA files
            with open(os.path.join(folder_path, filename), "r") as file:
                lines = file.readlines()

                # Separate the genes (assuming each gene is in its own header section)
                gene_stop_codons = 0
                current_gene_sequence = []

                for line in lines:
                    line = line.strip()
                    if line.startswith(">"):  # Header line, end of current gene
                        if current_gene_sequence:
                            # Count internal stop codons for the current gene
                            gene_stop_codons += count_internal_stop_codons("".join(current_gene_sequence))
                            current_gene_sequence = []  # Reset for the next gene
                    else:
                        current_gene_sequence.append(line)

                # Don't forget to count the last gene after the last header
                if current_gene_sequence:
                    gene_stop_codons += count_internal_stop_codons("".join(current_gene_sequence))

                # Append the result for this genome
                genome_stop_codons.append([filename, gene_stop_codons])

    return genome_stop_codons

# Function to save the results to a CSV file
def save_to_csv(genome_stop_codons, output_file):
    # Write the results to a CSV file
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Genome", "Internal Stop Codons"])
        writer.writerows(genome_stop_codons)

# Main function to execute the process
def main():
    # Set up argparse to handle input and output paths with flags
    parser = argparse.ArgumentParser(description="Count internal stop codons in gene sequences.")
    parser.add_argument("-i", "--input_folder", required=True, help="Path to the folder containing gene files (in FASTA format)")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output CSV file where results will be saved")

    # Parse the arguments
    args = parser.parse_args()

    # Process the folder and get the results
    genome_stop_codons = process_genome_folder(args.input_folder)

    # Save results to CSV
    save_to_csv(genome_stop_codons, args.output_file)
    print(f"Results saved to {args.output_file}")

if __name__ == "__main__":
    main()

