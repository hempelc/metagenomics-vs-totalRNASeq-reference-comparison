#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 22 Jul 2022

# This script processes coverage data from multiple replicates in .csv formats

# Note: only process either metagenomics or total rna-seq data at a time and adapt
# outfilename accordingly!

import glob
import os
import pandas as pd

# Full path to directory that contains .csv files of all samples to process
workdir = "/Users/christopherhempel/Google Drive/PhD UoG/Shea project/mapping_output_files/DNA"
outdir = "/Users/christopherhempel/Google Drive/PhD UoG/Shea project/data_processing_outdir"
# Prefix of output file
outfilename = "metagenomics"

sample_dfs_genome = {}
sample_dfs_ssu = {}
# Loop over all DNA/RNA files and save all as separate df
sample_csvs = glob.glob(os.path.join(workdir, "*.csv*"))
for file in sample_csvs:

    ## Read in file
    df = pd.read_csv(file).fillna(0)

    # Replace _ with - and one species name and make sure all the right tools are separated
    df["combo"] = df["combo"].str.replace("trimmed_at_phred_", "").str.replace(".fa", "").str.replace("IDBA_", "IDBA-").str.lower()
    df["species"] = df["species"].str.replace("_", " ").str.replace("Lactobacillus", "Limosilactobacillus")

    ## Separate df
    df_genome = df[df["level"]=="Genomes"]
    df_ssu = df[df["level"]=="ssuRNA_consensus"]

    ## Turn into heatmap format
    df_genome = df_genome.pivot_table(index="combo", columns="species", values='coverage[%]').reset_index()
    df_ssu = df_ssu.pivot_table(index="combo", columns="species", values='coverage[%]').reset_index()

    ## Generate tool columns and tidy up
    df_genome[["trimmingscore", "rrnasortingtool", "assemblytool"]] = df_genome["combo"].str.split("_", n=-1, expand=True)
    df_ssu[["trimmingscore", "rrnasortingtool", "assemblytool"]] = df_ssu["combo"].str.split("_", n=-1, expand=True)
    df_genome.columns.name = None
    df_ssu.columns.name = None
    df_genome = df_genome.drop(["combo"], axis=1)
    df_ssu = df_ssu.drop(["combo"], axis=1)

    ## Sort values
    df_genome = df_genome.sort_values(["assemblytool", "rrnasortingtool", "trimmingscore"])
    df_ssu = df_ssu.sort_values(["assemblytool", "rrnasortingtool", "trimmingscore"])

    ## Keep tool info and drop them so that dfs consist of only numbers
    tools = df_genome[["assemblytool", "rrnasortingtool", "trimmingscore"]]
    df_genome = df_genome.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1)
    df_ssu = df_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1)

    # Sort columns from most abundant to least abundant
    df_genome = df_genome.reindex(columns=['Listeria monocytogenes','Pseudomonas aeruginosa',
        'Bacillus subtilis','Saccharomyces cerevisiae','Salmonella enterica',
        'Escherichia coli','Limosilactobacillus fermentum','Enterococcus faecalis',
        'Cryptococcus neoformans','Staphylococcus aureus'])
    df_ssu = df_ssu.reindex(columns=['Listeria monocytogenes','Pseudomonas aeruginosa',
        'Bacillus subtilis','Saccharomyces cerevisiae','Salmonella enterica',
        'Escherichia coli','Limosilactobacillus fermentum','Enterococcus faecalis',
        'Cryptococcus neoformans','Staphylococcus aureus'])

    ## Save
    sample_dfs_genome[file] = df_genome
    sample_dfs_ssu[file] = df_ssu

# Take average across samples
df_average_genome = sample_dfs_genome[file] * 0
for sample_df_genome in sample_dfs_genome.values():
    df_average_genome = df_average_genome + sample_df_genome
df_average_genome = df_average_genome/len(sample_dfs_genome.values())

df_average_ssu = sample_dfs_ssu[file] * 0
for sample_df_genome in sample_dfs_ssu.values():
    df_average_ssu = df_average_ssu + sample_df_genome
df_average_ssu = df_average_ssu/len(sample_dfs_genome.values())

# Add tools
master_df_genome = pd.concat([df_average_genome, tools], axis=1)
master_df_ssu = pd.concat([df_average_ssu, tools], axis=1)

# Save
if not os.path.exists(outdir):
    os.makedirs(outdir)
master_df_genome.to_csv(os.path.join(outdir, outfilename + "_genome.csv"), index=False)
master_df_ssu.to_csv(os.path.join(outdir, outfilename + "_ssu.csv"), index=False)
