import os
import subprocess
import tempfile
import shutil
from Bio import SeqIO
import argparse
import csv

def count_sequences(fasta_file):
    return sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

def run_mmseqs_clustering(input_fasta, tmp_dir, output_dir, identity=0.90):
    db_path = os.path.join(tmp_dir, "seqDB")
    cluster_path = os.path.join(tmp_dir, "clusterRes")
    tsv_path = os.path.join(output_dir, "rep_seqs.tsv")

    subprocess.run(["mmseqs", "createdb", input_fasta, db_path], check=True)
    subprocess.run(["mmseqs", "cluster", db_path, cluster_path, tmp_dir, 
                    "-c", str(identity), "--min-seq-id", str(identity)], check=True)
    subprocess.run(["mmseqs", "createtsv", db_path, db_path, cluster_path, tsv_path], check=True)

    return tsv_path

def count_clusters(tsv_file):
    with open(tsv_file, "r") as f:
        clusters = set(line.split()[0] for line in f if line.strip())
    return len(clusters)

def calculate_redundancy_mmseqs(fasta_path, scratch_path):
    with tempfile.TemporaryDirectory(dir=scratch_path) as tmpdir:
        outdir = tempfile.mkdtemp(dir=scratch_path)
        total_seqs = count_sequences(fasta_path)
        tsv_file = run_mmseqs_clustering(fasta_path, tmpdir, outdir)
        num_clusters = count_clusters(tsv_file)
        shutil.rmtree(outdir)
        redundancy = 1 - (num_clusters / total_seqs) if total_seqs > 0 else 0
        return total_seqs, num_clusters, redundancy

def main(folder, scratch_path, output_csv):
    results = []

    for filename in os.listdir(folder):
        if filename.endswith(".faa"):
            filepath = os.path.join(folder, filename)
            try:
                total, clusters, red = calculate_redundancy_mmseqs(filepath, scratch_path)
                results.append((filename, total, clusters, round(red, 3)))
            except Exception as e:
                print(f"Error processing {filename}: {e}")

    # Save to CSV
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Filename", "Total_Sequences", "Clusters", "Redundancy"])
        writer.writerows(results)

    print(f"✅ Redundancy results saved to: {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate redundancy using MMseqs2 clustering at 95% identity.")
    parser.add_argument("folder", help="Folder containing FASTA files")
    parser.add_argument("--scratch", default="/tmp", help="Scratch directory for MMseqs temp files (default: /tmp)")
    parser.add_argument("--output", default="redundancy_results.csv", help="Output CSV filename")
    args = parser.parse_args()

    main(args.folder, args.scratch, args.output)

