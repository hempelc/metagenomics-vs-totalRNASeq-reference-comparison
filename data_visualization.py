import pandas as pd
import plotly.express as px
import statsmodels.api as sm

# Full path to specific .cvs files
metagenomics_genome_csv = "/Users/christopherhempel/Desktop/Pilot project results backup/pipeline_results_coverage"
metagenomics_ssu_csv = "/Users/christopherhempel/Desktop/Pilot project results backup/pipeline_results_coverage"
totalrnaseq_genome_csv = "/Users/christopherhempel/Desktop/Pilot project results backup/pipeline_results_coverage"
totalrnaseq_ssu_csv = "/Users/christopherhempel/Desktop/Pilot project results backup/pipeline_results_coverage"

# Species list
species = ['Listeria monocytogenes - 94.8%','Pseudomonas aeruginosa - 4.2%',
    'Bacillus subtilis - 0.7%','Saccharomyces cerevisiae - 0.23%','Escherichia coli - 0.089%',
    'Salmonella enterica - 0.059%','Lactobacillus fermentum - 0.015%','Enterococcus faecalis - 0.001%',
    'Cryptococcus neoformans - 0.00015%','Staphylococcus aureus - 0.0001%']

# Read in data
df_metagenomics_genome = pd.read_csv(metagenomics_genome_csv)
df_metagenomics_ssu = pd.read_csv(metagenomics_ssu_csv)
df_totalrnaseq_genome = pd.read_csv(totalrnaseq_genome_csv)
df_totalrnaseq_ssu = pd.read_csv(totalrnaseq_ssu_csv)

# Heatmap generation
heatmap_metagenomics_genome = px.imshow(df_metagenomics_genome.drop(["assemblytool", "rrnasortingtool", "trimmingscore"]),
    color_continuous_scale='tealgrn', labels=dict(x="Species (% Abundance)", y="Pipeline number", color="Genome coverage [%]"))
heatmap_metagenomics_ssu = px.imshow(df_metagenomics_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"]),
    color_continuous_scale='tealgrn', labels=dict(x="Species (% Abundance)", y="Pipeline number", color="SSU coverage [%]"))
heatmap_totalrnaseq_genome = px.imshow(df_totalrnaseq_genome.drop(["assemblytool", "rrnasortingtool", "trimmingscore"]),
    color_continuous_scale='tealgrn', labels=dict(x="Species (% Abundance)", y="Pipeline number", color="Genome coverage [%]"))
heatmap_totalrnaseq_ssu = px.imshow(df_totalrnaseq_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"]),
    color_continuous_scale='tealgrn', labels=dict(x="Species (% Abundance)", y="Pipeline number", color="Genome coverage [%]"))

heatmap_metagenomics_genome.show()
heatmap_metagenomics_ssu.show()
heatmap_totalrnaseq_genome.show()
heatmap_totalrnaseq_ssu.show()


# Correlation testing
## Add column for seqtype to test for correlation
df_metagenomics_genome["seqtype"] = ["metagenomics"] * len(df_metagenomics_genome.index)
df_metagenomics_ssu["seqtype"] = ["metagenomics"] * len(df_metagenomics_ssu.index)
df_totalrnaseq_genome["seqtype"] = ["totalrnaseq"] * len(df_totalrnaseq_genome.index)
df_totalrnaseq_ssu["seqtype"] = ["totalrnaseq"] * len(df_totalrnaseq_ssu.index)

## Conact dfs
df_genome = pd.concat([df_metagenomics_genome, df_totalrnaseq_genome])
df_ssu = pd.concat([df_metagenomics_ssu, df_totalrnaseq_ssu])

## Loop over species AND do all together
## STOPPED HERE
## Loop over species:
for spe in species:
### Set X and Y
X_genome = df_genome[["assemblytool", "rrnasortingtool", "trimmingscore", "seqtype"]]
Y_genome = df_genome[spe]
X_ssu = df_ssu[["assemblytool", "rrnasortingtool", "trimmingscore", "seqtype"]]
Y_ssu = df_ssu[spe]
X_genome_dummies = pd.get_dummies(X_genome)
X_genome_dummies = sm.add_constant(X_genome_dummies) # adding a constant
X_ssu_dummies = pd.get_dummies(X_ssu)
X_ssu_dummies = sm.add_constant(X_ssu_dummies) # adding a constant

# Train model
lr_model_genome = sm.OLS(Y_genome, X_genome_dummies).fit()
lr_model_ssu = sm.OLS(Y_genome, X_genome_dummies).fit()

## Summarize the output and extract coefs and p vals
lr_summary_genome = lr_model_genome.summary2().tables[1][['Coef.', 'P>|t|']]
lr_summary_genome = lr_summary_genome.rename(columns={"Coef.": "Coefficient", "P>|t|": "p-value"})
lr_summary_genome = lr_summary_genome.drop("const")
lr_summary_genome = lr_summary_genome.reset_index()
lr_summary_genome[['category', 'method']] = lr_summary_genome['index'].str.split('_', expand=True)
lr_summary_genome = lr_summary_genome.set_index('index')

lr_summary_genome.loc[lr_summary_genome["p-value"] <= 0.001 , 'significance_cat'] = "***"
lr_summary_genome.loc[lr_summary_genome["p-value"] > 0.001 , 'significance_cat'] = "**"
lr_summary_genome.loc[lr_summary_genome["p-value"] > 0.01 , 'significance_cat'] = "*"
lr_summary_genome.loc[lr_summary_genome["p-value"] > 0.05, 'significance_cat'] = ""

lr_summary_ssu = lr_model_ssu.summary2().tables[1][['Coef.', 'P>|t|']]
lr_summary_ssu = lr_summary_ssu.rename(columns={"Coef.": "Coefficient", "P>|t|": "p-value"})
lr_summary_ssu = lr_summary_ssu.drop("const")
lr_summary_ssu = lr_summary_ssu.reset_index()
lr_summary_ssu[['category', 'method']] = lr_summary_ssu['index'].str.split('_', expand=True)
lr_summary_ssu = lr_summary_ssu.set_index('index')

lr_summary_ssu.loc[lr_summary_ssu["p-value"] <= 0.001 , 'significance_cat'] = "***"
lr_summary_ssu.loc[lr_summary_ssu["p-value"] > 0.001 , 'significance_cat'] = "**"
lr_summary_ssu.loc[lr_summary_ssu["p-value"] > 0.01 , 'significance_cat'] = "*"
lr_summary_ssu.loc[lr_summary_ssu["p-value"] > 0.05, 'significance_cat'] = ""


# lr_summary_genome = lr_summary_genome.reindex(['rank_phylum', 'rank_class', 'rank_order', 'rank_family', 'rank_genus',
#     'rank_species', 'datatype_abundance', 'datatype_pa',
#     'seqtype_metagenomics', 'seqtype_totalrnaseq', 'model_knn',
#     'model_lor-lasso', 'model_lor-ridge', 'model_lsvc', 'model_mlp',
#     'model_rf', 'model_svc', 'model_xgb'])
## Plot
cols = ["#2995ff", "#d96c1e", "#4db444","#9f0281"]
bar_genome = px.bar(lr_summary_genome, x='method', y='Coefficient', color='category', text='significance_cat', color_discrete_sequence = cols)
bar_genome.update_layout(xaxis_tickangle=45)
bar_genome.update_traces(textposition='outside')
bar_genome.update_yaxes(nticks = 9)
bar_genome.show()
