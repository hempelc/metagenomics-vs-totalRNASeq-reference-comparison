import pandas as pd
import glob
import os

## ONLY ALL METAGENOMICS IN WORKDIR OR TOTALRNASEQ IN THE OTHER, ADAPT OUTFILENAME
# Full path to directory that contains .csv files of all samples to process
workdir = "/Users/christopherhempel/Google Drive/PhD UoG/Shea project/mapping_output_files/RNA"
outdir = "/Users/christopherhempel/Google Drive/PhD UoG/Shea project/data_processing_outdir"
# Prefix of output file
outfilename = "totalrnaseq"

# Define abundances
abundances = {'Bacillus subtilis': 0.7,
                    'Cryptococcus neoformans': 0.00015,
                    'Enterococcus faecalis': 0.001,
                    'Escherichia coli': 0.058,
                    'Limosilactobacillus fermentum': 0.015,
                    'Listeria monocytogenes': 94.8,
                    'Pseudomonas aeruginosa': 4.2,
                    'Saccharomyces cerevisiae': 0.23,
                    'Salmonella enterica': 0.059,
                    'Staphylococcus aureus': 0.0001}

# Loop over all DNA/RNA files and save all as separate df
sample_dfs_genome = {}
sample_dfs_ssu = {}

sample_csvs = glob.glob(os.path.join(workdir, "*.csv*"))
for file in sample_csvs:

    # Read in file
    df = pd.read_csv(file).fillna(0)

    # Replace _ with - and one species name and make sure all the right tools are separated
    df["combo"] = df["combo"].str.replace("trimmed_at_phred_", "").str.replace(".fa", "").str.replace("IDBA_", "IDBA-").str.lower()
    df["species"] = df["species"].str.replace("_", " ").str.replace("Lactobacillus", "Limosilactobacillus")

    # Add abundances to species for columns later
    for species, abundance in abundances.items():
        df["species"] = df["species"].str.replace(species, species + " - " + str(abundance) + "%")

    # Separate df
    df_genome = df[df["level"]=="Genomes"]
    df_ssu = df[df["level"]=="ssuRNA_consensus"]

    # Turn into heatmap format
    df_genome = df_genome.pivot_table(index="combo", columns="species", values='coverage[%]').reset_index()
    df_ssu = df_ssu.pivot_table(index="combo", columns="species", values='coverage[%]').reset_index()

    # Generate tool columns and tidy up
    df_genome[["trimmingscore", "rrnasortingtool", "assemblytool"]] = df_genome["combo"].str.split("_", n=-1, expand=True)
    df_ssu[["trimmingscore", "rrnasortingtool", "assemblytool"]] = df_ssu["combo"].str.split("_", n=-1, expand=True)
    df_genome.columns.name = None
    df_ssu.columns.name = None
    df_genome = df_genome.drop(["combo"], axis=1)
    df_ssu = df_ssu.drop(["combo"], axis=1)

    # Sort values
    df_genome = df_genome.sort_values(["assemblytool", "rrnasortingtool", "trimmingscore"])
    df_ssu = df_ssu.sort_values(["assemblytool", "rrnasortingtool", "trimmingscore"])

    # Keep tool info and drop them so that dfs consist of only numbers
    tools = df_genome[["assemblytool", "rrnasortingtool", "trimmingscore"]]
    df_genome = df_genome.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1)
    df_ssu = df_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1)

    # Sort columns from most abundant to least abundant
    df_genome = df_genome.reindex(columns=['Listeria monocytogenes - 94.8%','Pseudomonas aeruginosa - 4.2%',
        'Bacillus subtilis - 0.7%','Saccharomyces cerevisiae - 0.23%','Salmonella enterica - 0.059%',
        'Escherichia coli - 0.058%','Limosilactobacillus fermentum - 0.015%','Enterococcus faecalis - 0.001%',
        'Cryptococcus neoformans - 0.00015%','Staphylococcus aureus - 0.0001%'])
    df_ssu = df_ssu.reindex(columns=['Listeria monocytogenes - 94.8%','Pseudomonas aeruginosa - 4.2%',
        'Bacillus subtilis - 0.7%','Saccharomyces cerevisiae - 0.23%','Salmonella enterica - 0.059%',
        'Escherichia coli - 0.058%','Limosilactobacillus fermentum - 0.015%','Enterococcus faecalis - 0.001%',
        'Cryptococcus neoformans - 0.00015%','Staphylococcus aureus - 0.0001%'])

    # Save
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
