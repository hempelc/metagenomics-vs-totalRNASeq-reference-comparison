import pandas as pd
import plotly.express as px
import statsmodels.api as sm
from scipy import stats
import os
from scipy.stats import pointbiserialr #v1.7.3


workdir = "/Users/christopherhempel/Google Drive/PhD UoG/Shea project/data_processing_outdir/"
outdir = workdir

# Full path to specific .cvs files
metagenomics_genome_csv = os.path.join(workdir, "metagenomics_genome.csv")
metagenomics_ssu_csv = os.path.join(workdir, "metagenomics_ssu.csv")
totalrnaseq_genome_csv = os.path.join(workdir, "totalrnaseq_genome.csv")
totalrnaseq_ssu_csv = os.path.join(workdir, "totalrnaseq_ssu.csv")

# Reference df
cols=["#062a95","#cf8406"]
abundance_df = pd.DataFrame({"genomesize": [2992342, 6792330, 4045677, 12071326, 4809318, 4875441, 1905333, 2845392, 19051922, 2730326, 1000000, 5000000, 10000000, 15000000, 20000000],
    "abundance": [94.8, 4.2, 0.7, 0.23, 0.059, 0.058, 0.015, 0.001, 0.00015, 0.0001, 20, 20, 20, 20, 20],
    "species": ['Listeria monocytogenes', 'Pseudomonas aeruginosa', 'Bacillus subtilis', 'Saccharomyces cerevisiae', 'Salmonella enterica',
    'Escherichia coli', 'Limosilactobacillus fermentum', 'Enterococcus faecalis', 'Cryptococcus neoformans', 'Staphylococcus aureus', 'Staphylococcus aureus',
    'Cryptococcus neoformans', 'Enterococcus faecalis', 'Limosilactobacillus fermentum', 'Escherichia coli'],
    "domain": 3*["Bacteria"]+["Eukaryota"]+4*["Bacteria"]+["Eukaryota"]+["Bacteria"]+5*["ref"]})
reference = px.scatter(abundance_df, x="abundance", y="species", size="genomesize", text="abundance", color="domain", log_x=True, color_discrete_sequence=cols)
reference.update_traces(textposition='top center')
reference.update_yaxes(categoryorder='array', categoryarray= list(reversed(['Listeria monocytogenes', 'Pseudomonas aeruginosa', 'Bacillus subtilis', 'Saccharomyces cerevisiae', 'Salmonella enterica',
'Escherichia coli', 'Limosilactobacillus fermentum', 'Enterococcus faecalis', 'Cryptococcus neoformans', 'Staphylococcus aureus'])))
reference.show()
reference.write_image(os.path.join(outdir, "reference_bubbleplot.svg"))

# Save SSU gene lengths (although we don't use them)
ssulen = {'Bacillus subtilis': 1558,
    'Cryptococcus neoformans': 1513,
    'Enterococcus faecalis': 1562,
    'Escherichia coli': 1542,
    'Limosilactobacillus fermentum': 1578,
    'Listeria monocytogenes': 1552,
    'Pseudomonas aeruginosa': 1526,
    'Saccharomyces cerevisiae': 1553,
    'Salmonella enterica': 1534,
    'Staphylococcus aureus': 1556}

# Read in data
df_metagenomics_genome = pd.read_csv(metagenomics_genome_csv, dtype = {'trimmingscore': str})
df_metagenomics_ssu = pd.read_csv(metagenomics_ssu_csv, dtype = {'trimmingscore': str})
df_totalrnaseq_genome = pd.read_csv(totalrnaseq_genome_csv, dtype = {'trimmingscore': str})
df_totalrnaseq_ssu = pd.read_csv(totalrnaseq_ssu_csv, dtype = {'trimmingscore': str})


# ANOVAs
df_metagenomics_anova = pd.DataFrame()
df_totalrnaseq_anova = pd.DataFrame()
## Loop over species:
for spe in df_metagenomics_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).columns:
    # Metagenomics
    metagenomics_pval_trimmingscore = stats.f_oneway(df_metagenomics_ssu[spe][df_metagenomics_ssu['trimmingscore']=='5'],
      df_metagenomics_ssu[spe][df_metagenomics_ssu['trimmingscore']=='10'],
      df_metagenomics_ssu[spe][df_metagenomics_ssu['trimmingscore']=='15'],
      df_metagenomics_ssu[spe][df_metagenomics_ssu['trimmingscore']=='20'])[1]

    metagenomics_pval_rrnasortingtool = stats.f_oneway(df_metagenomics_ssu[spe][df_metagenomics_ssu['rrnasortingtool']=='barrnap'],
      df_metagenomics_ssu[spe][df_metagenomics_ssu['rrnasortingtool']=='rrnafilter'],
      df_metagenomics_ssu[spe][df_metagenomics_ssu['rrnasortingtool']=='sortmerna'],
      df_metagenomics_ssu[spe][df_metagenomics_ssu['rrnasortingtool']=='unsorted'])[1]

    metagenomics_pval_assemblytool = stats.f_oneway(df_metagenomics_ssu[spe][df_metagenomics_ssu['assemblytool']=='idba-tran'],
      df_metagenomics_ssu[spe][df_metagenomics_ssu['assemblytool']=='idba-ud'],
      df_metagenomics_ssu[spe][df_metagenomics_ssu['assemblytool']=='megahit'],
      df_metagenomics_ssu[spe][df_metagenomics_ssu['assemblytool']=='metaspades'],
      df_metagenomics_ssu[spe][df_metagenomics_ssu['assemblytool']=='rnaspades'],
      df_metagenomics_ssu[spe][df_metagenomics_ssu['assemblytool']=='spades'],
      df_metagenomics_ssu[spe][df_metagenomics_ssu['assemblytool']=='transabyss'])[1]

    # Total RNA Seq
    totalrnaseq_pval_trimmingscore = stats.f_oneway(df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['trimmingscore']=='5'],
      df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['trimmingscore']=='10'],
      df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['trimmingscore']=='15'],
      df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['trimmingscore']=='20'])[1]

    totalrnaseq_pval_rrnasortingtool = stats.f_oneway(df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['rrnasortingtool']=='barrnap'],
      df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['rrnasortingtool']=='rrnafilter'],
      df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['rrnasortingtool']=='sortmerna'],
      df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['rrnasortingtool']=='unsorted'])[1]

    totalrnaseq_pval_assemblytool = stats.f_oneway(df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['assemblytool']=='idba-tran'],
      df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['assemblytool']=='idba-ud'],
      df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['assemblytool']=='megahit'],
      df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['assemblytool']=='metaspades'],
      df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['assemblytool']=='rnaspades'],
      df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['assemblytool']=='spades'],
      df_totalrnaseq_ssu[spe][df_totalrnaseq_ssu['assemblytool']=='transabyss'])[1]

    df_metagenomics_anova[spe] = [metagenomics_pval_trimmingscore, metagenomics_pval_rrnasortingtool, metagenomics_pval_assemblytool]
    df_totalrnaseq_anova[spe] = [totalrnaseq_pval_trimmingscore, totalrnaseq_pval_rrnasortingtool, totalrnaseq_pval_assemblytool]
    df_metagenomics_anova.index = ["Trimming & quality filtering", "rRNA sorting", "Assembly"]
    df_totalrnaseq_anova.index = ["Trimming PHRED score", "rRNA sorting tool", "Assembler"]

## Plot
heatmap_metagenomics_anova = px.imshow(df_metagenomics_anova.transpose(), labels=dict(x="Step", y="Mock community species", color="P-value"), zmin=0, zmax=1, text_auto=".2f", color_continuous_scale="Blues_r")
heatmap_metagenomics_anova.update_xaxes(tickangle=35)
heatmap_metagenomics_anova.show()
heatmap_metagenomics_anova.write_image(os.path.join(outdir, "anova_metagenomics.svg"), height=650, width=600)
heatmap_metagenomics_anova.write_image(os.path.join(outdir, "anova_metagenomics.png"), height=650, width=600)

heatmap_totalrnaseq_anova = px.imshow(df_totalrnaseq_anova.transpose(), labels=dict(x="Step", y="Mock community species", color="P-value"), zmin=0, zmax=1, text_auto=".2f", color_continuous_scale="Blues_r")
heatmap_totalrnaseq_anova.update_xaxes(tickangle=35)
heatmap_totalrnaseq_anova.show()
heatmap_totalrnaseq_anova.write_image(os.path.join(outdir, "anova_totalrnaseq.svg"), height=650, width=600)
heatmap_totalrnaseq_anova.write_image(os.path.join(outdir, "anova_totalrnaseq.png"), height=650, width=600)


# Drop trimming score
df_metagenomics_genome = df_metagenomics_genome[df_metagenomics_genome['trimmingscore']=='20']
df_totalrnaseq_genome = df_totalrnaseq_genome[df_totalrnaseq_genome['trimmingscore']=='20']
df_metagenomics_ssu = df_metagenomics_ssu[df_metagenomics_ssu['trimmingscore']=='20']
df_totalrnaseq_ssu = df_totalrnaseq_ssu[df_totalrnaseq_ssu['trimmingscore']=='20']


# Correlation testing
df_metagenomics_cor = pd.DataFrame()
df_totalrnaseq_cor = pd.DataFrame()
## Loop over species:
for spe in df_metagenomics_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).columns:
    ### Set X and Y
    X_metagenomics = df_metagenomics_ssu[["assemblytool", "rrnasortingtool"]]
    X_totalrnaseq = df_totalrnaseq_ssu[["assemblytool", "rrnasortingtool"]]
    Y_metagenomics = df_metagenomics_ssu[spe]
    Y_totalrnaseq = df_totalrnaseq_ssu[spe]
    X_metagenomics_dummies = pd.get_dummies(X_metagenomics)
    X_metagenomics_dummies = sm.add_constant(X_metagenomics_dummies) # adding a constant
    X_totalrnaseq_dummies = pd.get_dummies(X_totalrnaseq)
    X_totalrnaseq_dummies = sm.add_constant(X_totalrnaseq_dummies) # adding a constant

    # Metagenomics
    cor_dic_metagenomics={}
    for col in X_metagenomics_dummies.columns:
        Y=Y_metagenomics.to_numpy()
        X=np.array(X_metagenomics_dummies[col])
        cor=pointbiserialr(X,Y)
        cor_dic_metagenomics[col]=cor
    cor_df_metagenomics = pd.DataFrame(cor_dic_metagenomics , index=["coefficient", "p-value"]).transpose().reset_index()
    cor_df_metagenomics = cor_df_metagenomics.drop(0)

    cor_df_metagenomics["species"] = [spe]*len(cor_df_metagenomics)

    cor_df_metagenomics.loc[cor_df_metagenomics["p-value"] <= 0.001 , 'significance_cat'] = 1
    cor_df_metagenomics.loc[cor_df_metagenomics["p-value"] > 0.001 , 'significance_cat'] = 0.7
    cor_df_metagenomics.loc[cor_df_metagenomics["p-value"] > 0.01 , 'significance_cat'] = 0.4
    cor_df_metagenomics.loc[cor_df_metagenomics["p-value"] > 0.05, 'significance_cat'] = 0.1

    df_metagenomics_cor = pd.concat([df_metagenomics_cor, cor_df_metagenomics])

    # Total RNA-Seq
    cor_dic_totalrnaseq={}
    for col in X_totalrnaseq_dummies.columns:
        Y=Y_totalrnaseq.to_numpy()
        X=np.array(X_totalrnaseq_dummies[col])
        cor=pointbiserialr(X,Y)
        cor_dic_totalrnaseq[col]=cor
    cor_df_totalrnaseq = pd.DataFrame(cor_dic_totalrnaseq , index=["coefficient", "p-value"]).transpose().reset_index()
    cor_df_totalrnaseq = cor_df_totalrnaseq.drop(0)

    cor_df_totalrnaseq["species"] = [spe]*len(cor_df_totalrnaseq)

    cor_df_totalrnaseq.loc[cor_df_totalrnaseq["p-value"] <= 0.001 , 'significance_cat'] = 1
    cor_df_totalrnaseq.loc[cor_df_totalrnaseq["p-value"] > 0.001 , 'significance_cat'] = 0.6
    cor_df_totalrnaseq.loc[cor_df_totalrnaseq["p-value"] > 0.01 , 'significance_cat'] = 0.35
    cor_df_totalrnaseq.loc[cor_df_totalrnaseq["p-value"] > 0.05, 'significance_cat'] = 0.1

    df_totalrnaseq_cor = pd.concat([df_totalrnaseq_cor, cor_df_totalrnaseq])


# Plot
bubble_metagenomcis_cor = px.scatter(df_metagenomics_cor, x="index", y="species", size="significance_cat",
    color="coefficient", height=550, width=700, range_color=(-1,1),
    color_continuous_scale="RdBu")
bubble_metagenomcis_cor.update_yaxes(categoryorder='array', categoryarray= list(reversed(['Listeria monocytogenes', 'Pseudomonas aeruginosa', 'Bacillus subtilis', 'Saccharomyces cerevisiae', 'Salmonella enterica',
'Escherichia coli', 'Limosilactobacillus fermentum', 'Enterococcus faecalis', 'Cryptococcus neoformans', 'Staphylococcus aureus'])))
bubble_metagenomcis_cor.update_xaxes(categoryorder='array', categoryarray= list(reversed(['assemblytool_idba-tran', 'assemblytool_idba-ud', 'assemblytool_megahit', 'assemblytool_metaspades', 'assemblytool_rnaspades',
'assemblytool_spades', 'assemblytool_transabyss', 'rrnasortingtool_barrnap', 'rrnasortingtool_rrnafilter', 'rrnasortingtool_sortmerna', 'rrnasortingtool_unsorted'])))
bubble_metagenomcis_cor.show()
bubble_metagenomcis_cor.write_image(os.path.join(outdir, "correlation_metagenomics.svg"))
bubble_metagenomcis_cor.write_image(os.path.join(outdir, "correlation_metagenomics.png"))

bubble_totalrnaseq_cor = px.scatter(df_totalrnaseq_cor, x="index", y="species", size="significance_cat",
    color="coefficient", height=550, width=700, range_color=(-1,1),
    color_continuous_scale="RdBu")
bubble_totalrnaseq_cor.update_yaxes(categoryorder='array', categoryarray= list(reversed(['Listeria monocytogenes', 'Pseudomonas aeruginosa', 'Bacillus subtilis', 'Saccharomyces cerevisiae', 'Salmonella enterica',
'Escherichia coli', 'Limosilactobacillus fermentum', 'Enterococcus faecalis', 'Cryptococcus neoformans', 'Staphylococcus aureus'])))
bubble_totalrnaseq_cor.update_xaxes(categoryorder='array', categoryarray= list(reversed(['assemblytool_idba-tran', 'assemblytool_idba-ud', 'assemblytool_megahit', 'assemblytool_metaspades', 'assemblytool_rnaspades',
'assemblytool_spades', 'assemblytool_transabyss', 'rrnasortingtool_barrnap', 'rrnasortingtool_rrnafilter', 'rrnasortingtool_sortmerna', 'rrnasortingtool_unsorted'])))
bubble_totalrnaseq_cor.show()
bubble_totalrnaseq_cor.write_image(os.path.join(outdir, "correlation_totalrnaseq.svg"))
bubble_totalrnaseq_cor.write_image(os.path.join(outdir, "correlation_totalrnaseq.png"))

# Coverage heatmaps
# Resort dfs
df_metagenomics_genome = df_metagenomics_genome.sort_values(["rrnasortingtool", "assemblytool"])
df_totalrnaseq_genome = df_totalrnaseq_genome.sort_values(["rrnasortingtool", "assemblytool"])
df_metagenomics_ssu = df_metagenomics_ssu.sort_values(["rrnasortingtool", "assemblytool"])

# Set index for readibility in plot
df_metagenomics_genome.index = df_metagenomics_genome["rrnasortingtool"] + "_" + df_metagenomics_genome["assemblytool"]
df_totalrnaseq_genome.index = df_totalrnaseq_genome["rrnasortingtool"] + "_" + df_totalrnaseq_genome["assemblytool"]
df_metagenomics_ssu.index = df_metagenomics_ssu["rrnasortingtool"] + "_" + df_metagenomics_ssu["assemblytool"]
df_totalrnaseq_ssu.index = df_totalrnaseq_ssu["assemblytool"] + "_" + df_totalrnaseq_ssu["rrnasortingtool"]

# Heatmap generation
heatmap_metagenomics_genome = px.imshow(df_metagenomics_genome.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).transpose(), aspect="auto", width=1150, labels=dict(x="Species (% Abundance)", y="Pipeline number", color="Genome coverage [%]"), zmin=0, zmax=100)
heatmap_metagenomics_ssu = px.imshow(df_metagenomics_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).transpose(), aspect="auto", width=1150, labels=dict(x="Species (% Abundance)", y="Pipeline number", color="SSU coverage [%]"), zmin=0, zmax=100)
heatmap_totalrnaseq_genome = px.imshow(df_totalrnaseq_genome.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).transpose(), aspect="auto", width=1150, labels=dict(x="Species (% Abundance)", y="Pipeline number", color="Genome coverage [%]"), zmin=0, zmax=100)
heatmap_totalrnaseq_ssu = px.imshow(df_totalrnaseq_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).transpose(), aspect="auto", width=1150, labels=dict(x="Species (% Abundance)", y="Pipeline number", color="SSU coverage [%]"), zmin=0, zmax=100)

# Plot
heatmap_metagenomics_genome.show()
heatmap_metagenomics_ssu.show()
heatmap_totalrnaseq_genome.show()
heatmap_totalrnaseq_ssu.show()

heatmap_metagenomics_genome.write_image(os.path.join(outdir, "coverage_metagenomics_genome.svg"))
heatmap_metagenomics_ssu.write_image(os.path.join(outdir, "coverage_metagenomics_ssu.svg"))
heatmap_totalrnaseq_genome.write_image(os.path.join(outdir, "coverage_totalrnaseq_genome.svg"))
heatmap_totalrnaseq_ssu.write_image(os.path.join(outdir, "coverage_totalrnaseq_ssu.svg"))

heatmap_metagenomics_genome.write_image(os.path.join(outdir, "coverage_metagenomics_genome.png"))
heatmap_metagenomics_ssu.write_image(os.path.join(outdir, "coverage_metagenomics_ssu.png"))
heatmap_totalrnaseq_genome.write_image(os.path.join(outdir, "coverage_totalrnaseq_genome.png"))
heatmap_totalrnaseq_ssu.write_image(os.path.join(outdir, "coverage_totalrnaseq_ssu.png"))
