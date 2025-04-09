import os
import subprocess
import argparse
from Bio import SeqIO
from math import ceil
import shutil
import tempfile

def split_fasta(input_path, chunk_dir, chunk_size=100):
    records = list(SeqIO.parse(input_path, "fasta"))
    total_chunks = ceil(len(records) / chunk_size)
    chunk_paths = []

    for i in range(total_chunks):
        chunk_records = records[i * chunk_size:(i + 1) * chunk_size]
        chunk_path = os.path.join(chunk_dir, f"chunk_{i:03d}.fasta")
        with open(chunk_path, "w") as out_fasta:
            SeqIO.write(chunk_records, out_fasta, "fasta")
        chunk_paths.append(chunk_path)

    return chunk_paths

def rename_genes(fasta_path, renamed_path, map_path):
    name_map = {}
    with open(renamed_path, "w") as out_fasta, open(map_path, "w") as map_file:
        for i, record in enumerate(SeqIO.parse(fasta_path, "fasta")):
            short_id = f"gene{i+1:04d}"
            name_map[short_id] = record.id
            map_file.write(f"{short_id}\t{record.id}\n")
            record.id = short_id
            record.description = ""
            SeqIO.write(record, out_fasta, "fasta")
    return name_map



def run_codonw(fasta_file, codonw_dir):
    # Short path to avoid buffer overflow
    short_fasta = os.path.join(codonw_dir, "temp.fasta")
    shutil.copy(fasta_file, short_fasta)

    subprocess.run(["codonw", "temp.fasta", "-enc", "-nomenu", "-silent"], cwd=codonw_dir, stdout=subprocess.DEVNULL)

    return os.path.join(codonw_dir, "temp.out")


def parse_enc_file(enc_file, name_map):
    enc_data = []
    with open(enc_file) as f_in:
        lines = f_in.readlines()
        for line in lines[1:]:  # Skip header
            parts = line.strip().split()
            if len(parts) >= 2:
                short_id = parts[0]
                enc_value = parts[1]
                original_id = name_map.get(short_id, short_id)
                enc_data.append((original_id, enc_value))
    return enc_data

def process_fasta(input_path, output_path, temp_dir, chunk_size=100):
    genome_id = os.path.splitext(os.path.basename(input_path))[0]
    genome_temp = os.path.join(temp_dir, genome_id)
    os.makedirs(genome_temp, exist_ok=True)

    chunk_paths = split_fasta(input_path, genome_temp, chunk_size)
    all_enc_data = []

    for i, chunk in enumerate(chunk_paths):
        renamed = os.path.join(genome_temp, f"chunk_{i:03d}_renamed.fasta")
        map_path = os.path.join(genome_temp, f"chunk_{i:03d}_map.tsv")

        name_map = rename_genes(chunk, renamed, map_path)

        # Use short CodonW path
        enc_file = run_codonw(renamed, genome_temp)
        enc_data = parse_enc_file(enc_file, name_map)
        all_enc_data.extend(enc_data)

        # Clean up
        os.remove(enc_file)
        os.remove(os.path.join(genome_temp, "temp.fasta"))

    # Write final merged output
    with open(output_path, "w") as f_out:
        for gene, enc in all_enc_data:
            f_out.write(f"{gene},{enc}\n")


def process_all(input_dir, output_dir, chunk_size=100):
    os.makedirs(output_dir, exist_ok=True)
    temp_dir = os.path.join(output_dir, "temp")
    os.makedirs(temp_dir, exist_ok=True)

    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta") or filename.endswith(".fna"):
            input_path = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, os.path.splitext(filename)[0] + ".out")
            print(f"Processing {filename}...")
            process_fasta(input_path, output_file, temp_dir, chunk_size)

    print("All genomes processed!")

def main():
    parser = argparse.ArgumentParser(description="Calculate ENC values using CodonW for multiple genomes with splitting.")
    parser.add_argument("-i", "--input_dir", required=True, help="Folder with input FASTA files (one per genome)")
    parser.add_argument("-o", "--output_dir", required=True, help="Folder to store output .out files")
    parser.add_argument("-c", "--chunk_size", type=int, default=100, help="Number of genes per CodonW run (default: 100)")
    args = parser.parse_args()

    process_all(args.input_dir, args.output_dir, args.chunk_size)

if __name__ == "__main__":
    main()

