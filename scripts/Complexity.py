import os
import math
import csv
import argparse

# Function to calculate Shannon entropy (sequence complexity)
def calculate_entropy(sequence):
    # Remove non-nucleotide characters (e.g., newlines, headers)
    sequence = ''.join([base for base in sequence if base in 'ATCG'])

    # Calculate frequencies of nucleotides
    nucleotide_count = {base: sequence.count(base) for base in 'ATCG'}
    sequence_length = len(sequence)

    # Calculate probabilities
    probabilities = [count / sequence_length for count in nucleotide_count.values()]

    # Calculate Shannon entropy
    entropy = -sum(p * math.log2(p) for p in probabilities if p > 0)
    return entropy

# Function to process each genome file in the folder
def process_genome_folder(folder_path):
    genome_complexity = []

    # Loop through all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".fasta") or filename.endswith(".fa")or filename.endswith(".fna"):  # Only process FASTA files
            with open(os.path.join(folder_path, filename), "r") as file:
                lines = file.readlines()

                # Join sequence lines and remove non-sequence characters (like header lines)
                sequence = "".join([line.strip() for line in lines if not line.startswith(">")])

                # Calculate sequence complexity (Shannon entropy)
                entropy = calculate_entropy(sequence)

                # Append the result to the list
                genome_complexity.append([filename, entropy])

    return genome_complexity

# Function to save the results to a CSV file
def save_to_csv(genome_complexity, output_file):
    # Write the results to a CSV file
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Genome", "Sequence Complexity (Entropy)"])
        writer.writerows(genome_complexity)

# Main function to execute the process
def main():
    # Set up argparse to handle input and output paths with flags
    parser = argparse.ArgumentParser(description="Calculate sequence complexity (Shannon entropy) for virus genomes.")
    parser.add_argument("-i", "--input_folder", required=True, help="Path to the folder containing genome files (in FASTA format)")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output CSV file where results will be saved")

    # Parse the arguments
    args = parser.parse_args()

    # Process the folder and get the results
    genome_complexity = process_genome_folder(args.input_folder)

    # Save results to CSV
    save_to_csv(genome_complexity, args.output_file)
    print(f"Results saved to {args.output_file}")

if __name__ == "__main__":
    main()

