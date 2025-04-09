
# GeneStat
![genestat (1)](https://github.com/user-attachments/assets/d2f22373-2210-422c-9266-d9835cfcb2a4)


GeneStat is a comprehensive Python pipeline that calculates various genome-based statistics, including:

- Basic genome assembly metrics (e.g., N50, GC content)
- Predicted gene count and average gene length using Prodigal
- Sequence complexity
- Codon usage statistics and principal components
- Effective Number of Codons (ENC) per gene and per genome using CodonW
- Internal stop codons
- Tetranucleotide frequency (TNF) and principal components

## 📦 Requirements

Before running the pipeline, make sure the following tools and libraries are installed:

### System Dependencies

You must have the following installed on your system and accessible from the command line:

- [Prodigal](https://github.com/hyattpd/Prodigal)
- [CodonW](http://codonw.sourceforge.net/)
- [SeqKit](https://bioinf.shenwei.me/seqkit/)

You can install them via:

```bash
# Prodigal
conda install -c bioconda prodigal

# CodonW
conda install -c bioconda codonw

# SeqKit
conda install -c bioconda seqkit
```

### Python Dependencies

Install the required Python libraries:

```bash
pip install biopython rich pandas scikit-learn
```

## 🧪 Usage

```
python genestat.py -i input_folder -o output_folder
```

- `input_folder`: Folder containing input FASTA files (one per genome)
- `output_folder`: Folder where output files and statistics will be saved

## 🔍 Pipeline Overview

1. **Prodigal ORF Prediction**: Predicts ORFs from each genome using a custom launcher.
2. **TNF Profiling**: Profiles tetranucleotide frequency for each genome.
3. **Complexity Calculation**: Measures entropy-based sequence complexity.
4. **Codon Usage**: Calculates codon usage frequency per genome.
5. **ENC Calculation**: Uses CodonW to estimate ENC values per gene and genome.
6. **Internal Stop Codons**: Detects internal stop codons in predicted ORFs.
7. **Basic Stats**: Uses SeqKit to summarize sequence assembly stats.
8. **Final Compilation**: All statistics are merged into a master CSV file called `All_genome_stats.csv`.

## 🗂 Output Files

- `ENC_per_gene.csv`: ENC value for each gene
- `ENC_genome_mean.csv`: Average ENC per genome
- `internal_stop_codons.csv`: Count of internal stop codons per genome
- `basic_stats.tsv`: SeqKit stats for original genomes
- `basic_stats_genes.tsv`: SeqKit stats for predicted genes
- `Complexity.csv`: Sequence entropy
- `codon_usage.csv`: Codon usage frequency
- `all_genomes_tetranucleotides_TNF_Profiles.csv`: TNF profiles
- `All_genome_stats.csv`: Combined metrics for each genome

## 🔧 Script Structure

This program assumes the following helper scripts are located in the `scripts/` folder:

- `prodigal_launcher.py`
- `TNF_profiler.py`
- `Complexity.py`
- `codon_usage.py`
- `calculate_enc.py`
- `internal_stop_codons.py`

Ensure they are executable and available in the same directory structure.

## 🧬 Example

```
python genestat.py -i example_data/ -o results/
```

This will generate all output files in the `results/` directory.

---

For any questions or contributions, feel free to open an issue or pull request.
