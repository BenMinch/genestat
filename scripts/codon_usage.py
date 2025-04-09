import os
import csv
import argparse
from collections import defaultdict
from Bio import SeqIO

CODONS = [a+b+c for a in "ATGC" for b in "ATGC" for c in "ATGC"]

def count_codons(seq):
    codon_count = defaultdict(int)
    seq = seq.upper().replace("\n", "").replace(" ", "")
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in CODONS:
            codon_count[codon] += 1
    return codon_count

def process_genome_file(filepath):
    codon_totals = defaultdict(int)
    total_codons = 0
    for record in SeqIO.parse(filepath, "fasta"):
        codon_count = count_codons(str(record.seq))
        for codon, count in codon_count.items():
            codon_totals[codon] += count
            total_codons += count
    codon_freq = {codon: (codon_totals[codon] / total_codons if total_codons > 0 else 0) for codon in CODONS}
    return codon_freq

def main(input_folder, output_csv):
    genome_results = []

    for filename in os.listdir(input_folder):
        if filename.endswith(".fasta") or filename.endswith(".fa") or filename.endswith(".fna"):
            filepath = os.path.join(input_folder, filename)
            codon_freq = process_genome_file(filepath)
            codon_freq["Genome"] = filename
            genome_results.append(codon_freq)

    with open(output_csv, "w", newline="") as csvfile:
        fieldnames = ["Genome"] + CODONS
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in genome_results:
            writer.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate codon usage frequency for each genome in a folder.")
    parser.add_argument("-i", "--input", required=True, help="Path to input folder containing genome FASTA files.")
    parser.add_argument("-o", "--output", required=True, help="Path to output CSV file.")

    args = parser.parse_args()
    main(args.input, args.output)

