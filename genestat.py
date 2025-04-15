import pandas as pd
import sys,argparse,os,re
from Bio import SeqIO
from math import ceil
import subprocess
from rich import print


def main():
    parser = argparse.ArgumentParser(description="Calculate ENC values using CodonW for multiple genomes with splitting.")
    parser.add_argument("-i", "--input_dir", required=True, help="Folder with input FASTA files (one per genome)")
    parser.add_argument("-o", "--output_dir", required=True, help="Folder to store output .out files")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    ##Run prodigal##
    print('[bold red]**Running Prodigal**[bold red]')
    os.mkdir(output_dir+'/Predicted_ORFs')
    prodigal= 'python scripts/prodigal_launcher.py '+input_dir+' '+output_dir+'/Predicted_ORFs'
    subprocess.call(prodigal, shell=True)

    ##Run TNF_profiler
    print('[bold orange]**Profiling TNF frequencies**[bold orange]')
    tnf='python scripts/TNF_profiler.py -i '+input_dir+' -o TNF_Profiles -b'
    subprocess.call(tnf, shell=True)
    move='mv all_genomes_tetranucleotides_TNF_Profiles.csv '+output_dir
    subprocess.call(move, shell=True)
    remove= 'rm -r TNF_Profiles'
    subprocess.call(remove, shell=True)

    ##Run Complexity
    print('[bold yellow]**Calculating Complexity**[bold yellow]')
    complexity='python scripts/Complexity.py -i '+input_dir+' -o '+output_dir+'/Complexity.csv'
    subprocess.call(complexity, shell=True)

    ##Run codon usage
    print('[bold green]**Calculating Codon Usage**[bold green]')
    gene_folder=output_dir+'/Predicted_ORFs'
    codon_usage='python scripts/codon_usage.py -i '+gene_folder+' -o '+output_dir+'/codon_usage.csv'
    subprocess.call(codon_usage, shell=True)
    ##Run Redundancy
    print('[bold magenta]**Calculating Redundancy**[bold magenta]')
    redundancy='python scripts/redundancy.py '+gene_folder+' --output '+output_dir+'/redundancy.csv'
    subprocess.call(redundancy, shell=True,stderr=open('redundancy.error', 'w'),stdout=open('redundancy.out', 'w'))

    ##Run ENC
    print('[bold blue]**Calculating ENC**[bold blue]')
    os.mkdir(output_dir+'/ENC')
    enc_folder=output_dir+'/ENC'
    enc='python scripts/calculate_enc.py -i '+gene_folder+' -o '+enc_folder
    subprocess.call(enc, shell=True,stderr=open('enc.error', 'w'),stdout=open('enc.out', 'w'))
    #Combine all the .out files in the ENC folder and add header for Genome and ENC
    enc_files = os.listdir(enc_folder)
    enc_files = [f for f in enc_files if f.endswith('.out')]
    enc_df = pd.DataFrame()
    for file in enc_files:
        genome = file.split('.genes')[0]
        df = pd.read_csv(enc_folder+'/'+file, sep=',', header=None)
        df.columns = ['Gene', 'ENC']
        df['Genome'] = genome
        enc_df = pd.concat([enc_df, df], ignore_index=True)
    #save the combined dataframe
    enc_df.to_csv(output_dir+'/ENC_per_gene.csv', index=False)
    #group by genome and get mean ENC
    #remove rows with ***** as ENC
    enc_df = enc_df[enc_df['ENC'] != '*****']
    enc_df['ENC'] = enc_df['ENC'].astype(float)

    enc_mean = enc_df.groupby('Genome')['ENC'].mean().reset_index()
    enc_mean.to_csv(output_dir+'/ENC_genome_mean.csv', index=False)

    ##Run internal stop codons
    print('[bold purple]**Calculating Internal Stop Codons**[bold purple]')
    stop='python scripts/internal_stop_codons.py -i '+gene_folder+' -o '+output_dir+'/internal_stop_codons.csv'
    subprocess.call(stop, shell=True)

    ##Run basic statistics
    print('[bold black]**Calculating Basic Statistics**[bold black]')
    stats='seqkit stats '+input_dir+'/* -a -T > '+output_dir+'/basic_stats.tsv'
    subprocess.call(stats, shell=True)
    ##Run stats on the gene folder to see how many ORFs were predicted
    stats='seqkit stats '+gene_folder+'/*.fna -a -T > '+output_dir+'/basic_stats_genes.tsv'
    subprocess.call(stats, shell=True)


    ##Combining everything together
    master_df = pd.read_csv(output_dir+'/basic_stats.tsv', sep='\t')
    #replace foldername in first column
    master_df['file'] = master_df['file'].str.replace(input_dir+'/', '')
    #only keep file, num_seqs, sum_len, N50, and GC
    master_df = master_df[['file', 'num_seqs', 'sum_len', 'N50', 'GC(%)']]
    #rename columns to Genome, Contigs, Length, N50, GC
    master_df.columns = ['Genome', 'Contigs', 'Length', 'N50', 'GC']
    
    #gene stats
    gene_stats = pd.read_csv(output_dir+'/basic_stats_genes.tsv', sep='\t')
    gene_stats['file'] = gene_stats['file'].str.replace(gene_folder+'/', '')
    #replace .genes.fna with ''
    gene_stats['file'] = gene_stats['file'].str.replace('.genes.fna', '')
    master_df['Num_Genes'] = master_df['Genome'].map(gene_stats.set_index('file')['num_seqs'])
    master_df['average_gene_length'] = master_df['Genome'].map(gene_stats.set_index('file')['avg_len'])

    #Get complexity
    complexity = pd.read_csv(output_dir+'/Complexity.csv')
    master_df['Complexity'] = master_df['Genome'].map(complexity.set_index('Genome')['Sequence Complexity (Entropy)'])

    #Get internal stop codons
    stop = pd.read_csv(output_dir+'/internal_stop_codons.csv')
    #replace .genes.fna with ''
    stop['Genome'] = stop['Genome'].str.replace('.genes.fna', '')
    master_df['Internal_Stop_Codons'] = master_df['Genome'].map(stop.set_index('Genome')['Internal Stop Codons'])

    #Calculate PC1 and PC2 for TNF
    tnf= pd.read_csv(output_dir+'/all_genomes_tetranucleotides_TNF_Profiles.csv')
    tnf2= tnf.drop(columns=['Genome'])
    #calculate PC1 and PC2
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(tnf2)
    #Make a dataframe with Genome and PC1 and PC2
    pca_df = pd.DataFrame(pca_result, columns=['TNF_PC1', 'TNF_PC2'])
    pca_df['Genome'] = tnf['Genome']
    #replace TNF in Genome with nothing
    pca_df['Genome'] = pca_df['Genome'].str.replace('TNF', '')
    #merge with master_df
    master_df = master_df.merge(pca_df, on='Genome', how='left')

    #Do the same thing with Codon usage
    codon_usage = pd.read_csv(output_dir+'/codon_usage.csv')
    #replace .genes.fna with ''
    codon_usage['Genome'] = codon_usage['Genome'].str.replace('.genes.fna', '')
    codon_usage2= codon_usage.drop(columns=['Genome'])
    #calculate PC1 and PC2
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(codon_usage2)
    #Make a dataframe with Genome and PC1 and PC2
    pca_df = pd.DataFrame(pca_result, columns=['Codon_PC1', 'Codon_PC2'])
    pca_df['Genome'] = codon_usage['Genome']
    master_df = master_df.merge(pca_df, on='Genome', how='left')

    #Do ENC
    enc = pd.read_csv(output_dir+'/ENC_genome_mean.csv')
    #replace .genes.fna with ''
    enc['Genome'] = enc['Genome'].str.replace('.genes.fna', '')
    #merge with master_df
    master_df = master_df.merge(enc, on='Genome', how='left')

    #Do Redundancy
    redundancy = pd.read_csv(output_dir+'/redundancy.csv')
    #replace .genes.fna with ''
    redundancy['Filename'] = redundancy['Filename'].str.replace('.faa', '')
    #merge with master_df
    master_df['Redundancy'] = master_df['Genome'].map(redundancy.set_index('Filename')['Redundancy'])
    #save the final dataframe
    master_df.to_csv(output_dir+'/All_genome_stats.csv', index=False)





if __name__ == "__main__":
    main()
