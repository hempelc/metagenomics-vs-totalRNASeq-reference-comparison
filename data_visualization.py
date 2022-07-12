import pandas as pd
import plotly.express as px
import statsmodels.api as sm

# Full path to specific .cvs files
metagenomics_genome_csv = "/Users/christopherhempel/Google Drive/PhD UoG/Shea project/data_processing_outdir/metagenomics_genome.csv"
metagenomics_ssu_csv = "/Users/christopherhempel/Google Drive/PhD UoG/Shea project/data_processing_outdir/metagenomics_ssu.csv"
totalrnaseq_genome_csv = "/Users/christopherhempel/Google Drive/PhD UoG/Shea project/data_processing_outdir/totalrnaseq_genome.csv"
totalrnaseq_ssu_csv = "/Users/christopherhempel/Google Drive/PhD UoG/Shea project/data_processing_outdir/totalrnaseq_ssu.csv"

# Species list
species = ['Listeria monocytogenes - 94.8%','Pseudomonas aeruginosa - 4.2%',
    'Bacillus subtilis - 0.7%','Saccharomyces cerevisiae - 0.23%','Salmonella enterica - 0.059%',
    'Escherichia coli - 0.058%','Limosilactobacillus fermentum - 0.015%','Enterococcus faecalis - 0.001%',
    'Cryptococcus neoformans - 0.00015%','Staphylococcus aureus - 0.0001%']

# Define ref lengths
genomelen = {'Bacillus subtilis - 0.7%': 4045677,
    'Cryptococcus neoformans - 0.00015%': 19051922,
    'Enterococcus faecalis - 0.001%': 2845392,
    'Escherichia coli - 0.058%': 4875441,
    'Limosilactobacillus fermentum - 0.015%': 1905333,
    'Listeria monocytogenes - 94.8%': 2992342,
    'Pseudomonas aeruginosa - 4.2%': 6792330,
    'Saccharomyces cerevisiae - 0.23%': 12071326,
    'Salmonella enterica - 0.059%': 4809318,
    'Staphylococcus aureus - 0.0001%': 2730326}

genomelen_short = {'Bacillus subtilis - 0.7%': "4 Mb",
    'Cryptococcus neoformans - 0.00015%': "19 Mb",
    'Enterococcus faecalis - 0.001%': "3 Mb",
    'Escherichia coli - 0.058%': "5 Mb",
    'Limosilactobacillus fermentum - 0.015%': "2 Mb",
    'Listeria monocytogenes - 94.8%': "3 Mb",
    'Pseudomonas aeruginosa - 4.2%': "7 Mb",
    'Saccharomyces cerevisiae - 0.23%': "12 Mb",
    'Salmonella enterica - 0.059%': "5 Mb",
    'Staphylococcus aureus - 0.0001%': "3 Mb"}

ssulen = {'Bacillus subtilis - 0.7%': 1558,
    'Cryptococcus neoformans - 0.00015%': 1513,
    'Enterococcus faecalis - 0.001%': 1562,
    'Escherichia coli - 0.058%': 1542,
    'Limosilactobacillus fermentum - 0.015%': 1578,
    'Listeria monocytogenes - 94.8%': 1552,
    'Pseudomonas aeruginosa - 4.2%': 1526,
    'Saccharomyces cerevisiae - 0.23%': 1553,
    'Salmonella enterica - 0.059%': 1534,
    'Staphylococcus aureus - 0.0001%': 1556}

# Read in data
df_metagenomics_genome = pd.read_csv(metagenomics_genome_csv, dtype = {'trimmingscore': str})
df_metagenomics_ssu = pd.read_csv(metagenomics_ssu_csv, dtype = {'trimmingscore': str})
df_totalrnaseq_genome = pd.read_csv(totalrnaseq_genome_csv, dtype = {'trimmingscore': str})
df_totalrnaseq_ssu = pd.read_csv(totalrnaseq_ssu_csv, dtype = {'trimmingscore': str})

cols_genome = df_metagenomics_genome.columns
for species, length in genomelen_short.items():
    cols_genome = cols_genome.str.replace(species, species + " - " + str(length))

cols_ssu = df_metagenomics_genome.columns
for species, length in ssulen.items():
    cols_ssu = cols_ssu.str.replace(species, species + " - " + str(length))

df_metagenomics_genome.columns = cols_genome
df_totalrnaseq_genome.columns = cols_genome
df_metagenomics_ssu.columns = cols_ssu
df_totalrnaseq_ssu.columns = cols_ssu

# Resort dfs
df_metagenomics_genome = df_metagenomics_genome.sort_values(["rrnasortingtool", "assemblytool", "trimmingscore"])
df_totalrnaseq_genome = df_totalrnaseq_genome.sort_values(["rrnasortingtool", "assemblytool", "trimmingscore"])
df_metagenomics_genome = df_metagenomics_genome.reset_index().drop('index', axis=1)
df_totalrnaseq_genome = df_totalrnaseq_genome.reset_index().drop('index', axis=1)
df_metagenomics_ssu = df_metagenomics_ssu.sort_values(["rrnasortingtool", "assemblytool", "trimmingscore"])
df_metagenomics_ssu = df_metagenomics_ssu.reset_index().drop('index', axis=1)


# Heatmap generation
heatmap_metagenomics_genome = px.imshow(df_metagenomics_genome.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1), labels=dict(x="Species (% Abundance)", y="Pipeline number", color="Genome coverage [%]"), zmin=0, zmax=100)
heatmap_metagenomics_ssu = px.imshow(df_metagenomics_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1),labels=dict(x="Species (% Abundance)", y="Pipeline number", color="SSU coverage [%]"), zmin=0, zmax=100)
heatmap_totalrnaseq_genome = px.imshow(df_totalrnaseq_genome.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1),labels=dict(x="Species (% Abundance)", y="Pipeline number", color="Genome coverage [%]"), zmin=0, zmax=100)
heatmap_totalrnaseq_ssu = px.imshow(df_totalrnaseq_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1), labels=dict(x="Species (% Abundance)", y="Pipeline number", color="SSU coverage [%]"), zmin=0, zmax=100)

# Plot
heatmap_metagenomics_genome.show()
heatmap_metagenomics_ssu.show()
heatmap_totalrnaseq_genome.show()
heatmap_totalrnaseq_ssu.show()


# Correlation testing
## Loop over species
df_metagenomics = pd.DataFrame()
df_totalrnaseq = pd.DataFrame()

## Loop over species:
for spe in df_metagenomics_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).columns:
    ### Set X and Y
    X_metagenomics = df_metagenomics_ssu[["assemblytool", "rrnasortingtool", "trimmingscore"]]
    Y_metagenomics = df_metagenomics_ssu[spe]
    X_totalrnaseq = df_totalrnaseq_ssu[["assemblytool", "rrnasortingtool", "trimmingscore"]]
    Y_totalrnaseq = df_totalrnaseq_ssu[spe]
    X_metagenomics_dummies = pd.get_dummies(X_metagenomics)
    X_metagenomics_dummies = sm.add_constant(X_metagenomics_dummies) # adding a constant
    X_totalrnaseq_dummies = pd.get_dummies(X_totalrnaseq)
    X_totalrnaseq_dummies = sm.add_constant(X_totalrnaseq_dummies) # adding a constant

    # Train model
    lr_model_metagenomics = sm.OLS(Y_metagenomics, X_metagenomics_dummies).fit()
    lr_model_totalrnaseq = sm.OLS(Y_totalrnaseq, X_totalrnaseq_dummies).fit()

    ## Summarize the output and extract coefs and p vals
    lr_summary_metagenomics = lr_model_metagenomics.summary2().tables[1][['Coef.', 'P>|t|']]
    lr_summary_metagenomics = lr_summary_metagenomics.rename(columns={"Coef.": "Coefficient", "P>|t|": "p-value"})
    lr_summary_metagenomics = lr_summary_metagenomics.drop("const")
    lr_summary_metagenomics = lr_summary_metagenomics.reset_index()
    lr_summary_metagenomics[['category', 'method']] = lr_summary_metagenomics['index'].str.split('_', expand=True)
    lr_summary_metagenomics = lr_summary_metagenomics.set_index('index')

    lr_summary_metagenomics.loc[lr_summary_metagenomics["p-value"] <= 0.001 , 'significance_cat'] = 1
    lr_summary_metagenomics.loc[lr_summary_metagenomics["p-value"] > 0.001 , 'significance_cat'] = 0.5
    lr_summary_metagenomics.loc[lr_summary_metagenomics["p-value"] > 0.01 , 'significance_cat'] = 0.25
    lr_summary_metagenomics.loc[lr_summary_metagenomics["p-value"] > 0.05, 'significance_cat'] = 0.1

    lr_summary_metagenomics["species"] = [spe]*len(lr_summary_metagenomics)

    lr_summary_totalrnaseq = lr_model_totalrnaseq.summary2().tables[1][['Coef.', 'P>|t|']]
    lr_summary_totalrnaseq = lr_summary_totalrnaseq.rename(columns={"Coef.": "Coefficient", "P>|t|": "p-value"})
    lr_summary_totalrnaseq = lr_summary_totalrnaseq.drop("const")
    lr_summary_totalrnaseq = lr_summary_totalrnaseq.reset_index()
    lr_summary_totalrnaseq[['category', 'method']] = lr_summary_totalrnaseq['index'].str.split('_', expand=True)
    lr_summary_totalrnaseq = lr_summary_totalrnaseq.set_index('index')

    lr_summary_totalrnaseq.loc[lr_summary_totalrnaseq["p-value"] <= 0.001 , 'significance_cat'] = 1
    lr_summary_totalrnaseq.loc[lr_summary_totalrnaseq["p-value"] > 0.001 , 'significance_cat'] = 0.5
    lr_summary_totalrnaseq.loc[lr_summary_totalrnaseq["p-value"] > 0.01 , 'significance_cat'] = 0.25
    lr_summary_totalrnaseq.loc[lr_summary_totalrnaseq["p-value"] > 0.05, 'significance_cat'] = 0.1

    lr_summary_totalrnaseq["species"] = [spe]*len(lr_summary_totalrnaseq)

    df_metagenomics = pd.concat([df_metagenomics, lr_summary_metagenomics])
    df_totalrnaseq = pd.concat([df_totalrnaseq, lr_summary_totalrnaseq])

# Plot
maximum = df_metagenomics["Coefficient"].max()
minimum = df_metagenomics["Coefficient"].min()
scale_range = max(abs(maximum), abs(minimum))
fig = px.scatter(df_metagenomics, x="species", y="method", size="significance_cat",
    color="Coefficient", height=1000, range_color=(-scale_range,scale_range),
    color_continuous_scale="RdBu")
fig.show()

maximum = df_totalrnaseq["Coefficient"].max()
minimum = df_totalrnaseq["Coefficient"].min()
scale_range = max(abs(maximum), abs(minimum))
fig = px.scatter(df_totalrnaseq, x="species", y="method", size="significance_cat",
    color="Coefficient", height=1000, range_color=(-scale_range,scale_range),
    color_continuous_scale="RdBu")
fig.show()
