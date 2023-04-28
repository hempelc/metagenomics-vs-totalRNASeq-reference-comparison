#!/usr/bin/env python3

# Written by Christopher Hempel (hempelc@uoguelph.ca) on 22 Jul 2022

# This script performs statistical tests on SSU rRNA and genome coverage from
# microbial mock community species and visualizes results

import os
import pandas as pd
import plotly.express as px
import numpy as np
import statsmodels.api as sm
import plotly.graph_objects as go
from scipy import stats
from statsmodels.formula.api import ols
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import QuantileTransformer


# Manual parameters
workdir = "/Users/christopherhempel/Google Drive/PhD UoG/Shea project/data_processing_outdir/"
outdir = workdir
# Full path to .cvs files generated in data_processing.py script
metagenomics_genome_csv = os.path.join(workdir, "metagenomics_genome.csv")
metagenomics_ssu_csv = os.path.join(workdir, "metagenomics_ssu.csv")
totalrnaseq_genome_csv = os.path.join(workdir, "totalrnaseq_genome.csv")
totalrnaseq_ssu_csv = os.path.join(workdir, "totalrnaseq_ssu.csv")
# Full path to masterlist files for aquarium samples
coverage_tsv = os.path.join(workdir, "masterlist_coverage.tsv")
taxa_tsv = os.path.join(workdir, "masterlist_taxa.tsv")

# Make a reference df based on abundances and genome sizes given by ZymoResearch and plot it
cols = ["#062a95", "#cf8406"]
reference_df = pd.DataFrame({"genomesize": [2992342, 6792330, 4045677, 12071326, 4809318, 4875441, 1905333, 2845392, 19051922, 2730326, 1000000, 5000000, 10000000, 15000000, 20000000],
                             "abundance": [94.8, 4.2, 0.7, 0.23, 0.059, 0.058, 0.015, 0.001, 0.00015, 0.0001, 20, 20, 20, 20, 20],
                             "species": ['Listeria monocytogenes', 'Pseudomonas aeruginosa', 'Bacillus subtilis', 'Saccharomyces cerevisiae', 'Salmonella enterica',
                                         'Escherichia coli', 'Limosilactobacillus fermentum', 'Enterococcus faecalis', 'Cryptococcus neoformans', 'Staphylococcus aureus', 'Staphylococcus aureus',
                                         'Cryptococcus neoformans', 'Enterococcus faecalis', 'Limosilactobacillus fermentum', 'Escherichia coli'],
                             "domain": 3*["Bacteria"]+["Eukaryota"]+4*["Bacteria"]+["Eukaryota"]+["Bacteria"]+5*["ref"]})
reference = px.scatter(reference_df, x="abundance", y="species", size="genomesize",
                       text="genomesize", color="domain", log_x=True, color_discrete_sequence=cols)
reference.update_traces(textposition='top center')
reference.update_yaxes(categoryorder='array', categoryarray=list(reversed(['Listeria monocytogenes', 'Pseudomonas aeruginosa', 'Bacillus subtilis', 'Saccharomyces cerevisiae', 'Salmonella enterica',
                                                                           'Escherichia coli', 'Limosilactobacillus fermentum', 'Enterococcus faecalis', 'Cryptococcus neoformans', 'Staphylococcus aureus'])))
reference.show()
reference.write_image(os.path.join(outdir, "reference_bubbleplot.svg"))


# Save SSU gene lengths (we actually don't use them for anything but we initially
# planned to, so we keep this information, doesn't hurt)
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
df_metagenomics_genome = pd.read_csv(
    metagenomics_genome_csv, dtype={'trimmingscore': str})
df_metagenomics_ssu = pd.read_csv(
    metagenomics_ssu_csv, dtype={'trimmingscore': str})
df_totalrnaseq_genome = pd.read_csv(
    totalrnaseq_genome_csv, dtype={'trimmingscore': str})
df_totalrnaseq_ssu = pd.read_csv(
    totalrnaseq_ssu_csv, dtype={'trimmingscore': str})
df_coverage = pd.read_csv(coverage_tsv, sep="\t")
df_taxa = pd.read_csv(taxa_tsv, sep="\t")
df_taxa = df_taxa.sort_values(by="sample")


# # ANOVAs - note: replaced by linear model
# df_metagenomics_anova = pd.DataFrame()
# df_totalrnaseq_anova = pd.DataFrame()
# ## Loop over species:
# for spe in df_metagenomics_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).columns:
#     # Metagenomics
#     twowayanovadf_metagenomics = df_metagenomics_ssu[[spe, 'trimmingscore', 'assemblytool', 'rrnasortingtool']].rename(columns={spe: spe.replace(" ", "_")})
#     data = twowayanovadf_metagenomics[spe.replace(" ", "_")].to_numpy().reshape(-1,1)
#     # quantile transform the raw data
#     quantile = QuantileTransformer(output_distribution='normal')
#     data_trans = quantile.fit_transform(data)
#     twowayanovadf_metagenomics[spe.replace(" ", "_")] = data_trans.reshape(len(data_trans))
#
#
#     ## Performing two-way ANOVA
#     model = ols(spe.replace(" ", "_") + ' ~ C(trimmingscore) + C(rrnasortingtool) + \
#         C(assemblytool) + C(trimmingscore):C(rrnasortingtool) + C(trimmingscore):\
#         C(assemblytool) + C(rrnasortingtool):C(assemblytool)',
#         data=twowayanovadf_metagenomics).fit()
#     twowayanovaresult_metagenomics = sm.stats.anova_lm(model, typ=3)
#     metagenomics_pval_trimmingscore = twowayanovaresult_metagenomics['PR(>F)'][1]
#     metagenomics_pval_rrnasortingtool = twowayanovaresult_metagenomics['PR(>F)'][2]
#     metagenomics_pval_assemblytool = twowayanovaresult_metagenomics['PR(>F)'][3]
#     metagenomics_pval_interaction1 = twowayanovaresult_metagenomics['PR(>F)'][4]
#     metagenomics_pval_interaction2 = twowayanovaresult_metagenomics['PR(>F)'][5]
#     metagenomics_pval_interaction3 = twowayanovaresult_metagenomics['PR(>F)'][6]
#
#
#     ## One-way ANOVA
#     # model = ols(spe.replace(" ", "_") + ' ~ C(trimmingscore)', data=df_metagenomics_ssu[[spe, 'trimmingscore', 'assemblytool', 'rrnasortingtool']].rename(columns={spe: spe.replace(" ", "_")})).fit()
#     # aov_table = sm.stats.anova_lm(model, typ=2)
#     # metagenomics_pval_trimmingscore = aov_table['PR(>F)'][0]
#     #
#     # model = ols(spe.replace(" ", "_") + ' ~ C(rrnasortingtool)', data=df_metagenomics_ssu[[spe, 'trimmingscore', 'assemblytool', 'rrnasortingtool']].rename(columns={spe: spe.replace(" ", "_")})).fit()
#     # aov_table = sm.stats.anova_lm(model, typ=2)
#     # metagenomics_pval_rrnasortingtool = aov_table['PR(>F)'][0]
#     #
#     # model = ols(spe.replace(" ", "_") + ' ~ C(assemblytool)', data=df_metagenomics_ssu[[spe, 'trimmingscore', 'assemblytool', 'rrnasortingtool']].rename(columns={spe: spe.replace(" ", "_")})).fit()
#     # aov_table = sm.stats.anova_lm(model, typ=2)
#     # metagenomics_pval_assemblytool = aov_table['PR(>F)'][0]
#
#     # Total RNA Seq
#     twowayanovadf_totalrnaseq = df_totalrnaseq_ssu[[spe, 'trimmingscore', 'assemblytool', 'rrnasortingtool']].rename(columns={spe: spe.replace(" ", "_")})
#     data = twowayanovadf_totalrnaseq[spe.replace(" ", "_")].to_numpy().reshape(-1,1)
#     # quantile transform the raw data
#     quantile = QuantileTransformer(output_distribution='normal')
#     data_trans = quantile.fit_transform(data)
#     twowayanovadf_totalrnaseq[spe.replace(" ", "_")] = data_trans.reshape(len(data_trans))
#     ## Performing two-way ANOVA
#     model = ols(spe.replace(" ", "_") + ' ~ C(trimmingscore) + C(rrnasortingtool) + \
#         C(assemblytool) + C(trimmingscore):C(rrnasortingtool) + C(trimmingscore):\
#         C(assemblytool) + C(rrnasortingtool):C(assemblytool)',
#         data=twowayanovadf_totalrnaseq).fit()
#     twowayanovaresult_totalrnaseq = sm.stats.anova_lm(model, typ=3)
#     totalrnaseq_pval_trimmingscore = twowayanovaresult_totalrnaseq['PR(>F)'][1]
#     totalrnaseq_pval_rrnasortingtool = twowayanovaresult_totalrnaseq['PR(>F)'][2]
#     totalrnaseq_pval_assemblytool = twowayanovaresult_totalrnaseq['PR(>F)'][3]
#     totalrnaseq_pval_interaction1 = twowayanovaresult_totalrnaseq['PR(>F)'][4]
#     totalrnaseq_pval_interaction2 = twowayanovaresult_totalrnaseq['PR(>F)'][5]
#     totalrnaseq_pval_interaction3 = twowayanovaresult_totalrnaseq['PR(>F)'][6]
#
#     ## One-way ANOVA
#     # model = ols(spe.replace(" ", "_") + ' ~ C(trimmingscore)', data=df_totalrnaseq_ssu[[spe, 'trimmingscore', 'assemblytool', 'rrnasortingtool']].rename(columns={spe: spe.replace(" ", "_")})).fit()
#     # aov_table = sm.stats.anova_lm(model, typ=2)
#     # totalrnaseq_pval_trimmingscore = aov_table['PR(>F)'][0]
#     #
#     # model = ols(spe.replace(" ", "_") + ' ~ C(rrnasortingtool)', data=df_totalrnaseq_ssu[[spe, 'trimmingscore', 'assemblytool', 'rrnasortingtool']].rename(columns={spe: spe.replace(" ", "_")})).fit()
#     # aov_table = sm.stats.anova_lm(model, typ=2)
#     # totalrnaseq_pval_rrnasortingtool = aov_table['PR(>F)'][0]
#     #
#     # model = ols(spe.replace(" ", "_") + ' ~ C(assemblytool)', data=df_totalrnaseq_ssu[[spe, 'trimmingscore', 'assemblytool', 'rrnasortingtool']].rename(columns={spe: spe.replace(" ", "_")})).fit()
#     # aov_table = sm.stats.anova_lm(model, typ=2)
#     # totalrnaseq_pval_assemblytool = aov_table['PR(>F)'][0]
#
#     # Summarize information one-way ANOVA
#     # df_metagenomics_anova[spe] = [metagenomics_pval_trimmingscore, metagenomics_pval_rrnasortingtool, metagenomics_pval_assemblytool]
#     # df_totalrnaseq_anova[spe] = [totalrnaseq_pval_trimmingscore, totalrnaseq_pval_rrnasortingtool, totalrnaseq_pval_assemblytool]
#     # df_metagenomics_anova.index = ["Trimming & quality filtering", "rRNA sorting", "Assembly"]
#     # df_totalrnaseq_anova.index = ["Trimming & quality filtering", "rRNA sorting ", "Assembly"]
#
#
#     # Summarize information two-way ANOVA
#     df_metagenomics_anova[spe] = [metagenomics_pval_trimmingscore, metagenomics_pval_rrnasortingtool, metagenomics_pval_assemblytool, metagenomics_pval_interaction1, metagenomics_pval_interaction2, metagenomics_pval_interaction3]
#     df_totalrnaseq_anova[spe] = [totalrnaseq_pval_trimmingscore, totalrnaseq_pval_rrnasortingtool, totalrnaseq_pval_assemblytool, totalrnaseq_pval_interaction1, totalrnaseq_pval_interaction2, totalrnaseq_pval_interaction3]
#     df_metagenomics_anova.index = ["Trimming & quality filtering", "rRNA sorting", "Assembly", "Trimming & quality filtering:rRNA sorting", "Trimming & quality filtering:Assembly", "rRNA sorting:Assembly"]
#     df_totalrnaseq_anova.index = ["Trimming & quality filtering", "rRNA sorting", "Assembly", "Trimming & quality filtering:rRNA sorting", "Trimming & quality filtering:Assembly", "rRNA sorting:Assembly"]
#
# ## Plot
# heatmap_metagenomics_anova = px.imshow(df_metagenomics_anova.transpose(), labels=dict(x="Step", y="Mock community species", color="P-value"), zmin=0, zmax=1, text_auto=".2f", color_continuous_scale="Blues_r")
# heatmap_metagenomics_anova.update_xaxes(tickangle=35)
# heatmap_metagenomics_anova.show()
# heatmap_metagenomics_anova.write_image(os.path.join(outdir, "anova_metagenomics.svg"), height=650, width=600)
# heatmap_metagenomics_anova.write_image(os.path.join(outdir, "anova_metagenomics.png"), height=650, width=600)
#
# heatmap_totalrnaseq_anova = px.imshow(df_totalrnaseq_anova.transpose(), labels=dict(x="Step", y="Mock community species", color="P-value"), zmin=0, zmax=1, text_auto=".2f", color_continuous_scale="Blues_r")
# heatmap_totalrnaseq_anova.update_xaxes(tickangle=35)
# heatmap_totalrnaseq_anova.show()
# heatmap_totalrnaseq_anova.write_image(os.path.join(outdir, "anova_totalrnaseq.svg"), height=650, width=600)
# heatmap_totalrnaseq_anova.write_image(os.path.join(outdir, "anova_totalrnaseq.png"), height=650, width=600)


# Ordinary Least Squares linear regression
df_metagenomics_ols = pd.DataFrame()
df_totalrnaseq_ols = pd.DataFrame()
quantile = QuantileTransformer(output_distribution='normal')

# Loop over species:
for spe in df_metagenomics_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).columns:
    # Metagenomics
    # Prepare df
    ols_df_metagenomics = df_metagenomics_ssu[[
        spe, 'trimmingscore', 'assemblytool', 'rrnasortingtool']].rename(columns={spe: spe.replace(" ", "_")})
    ols_df_metagenomics_Y = ols_df_metagenomics[spe.replace(" ", "_")]
    # Set up Y. Y data is not normal distributed, so we do a quantile transformation to transform into a normal distribution
    ols_df_metagenomics_Y = ols_df_metagenomics_Y.to_numpy().reshape(-1, 1)
    ols_df_metagenomics_Y = quantile.fit_transform(ols_df_metagenomics_Y)
    ols_df_metagenomics_Y = ols_df_metagenomics_Y.reshape(
        len(ols_df_metagenomics_Y))
    # Set up X as dummies
    ols_df_metagenomics_X_dummies = pd.get_dummies(
        ols_df_metagenomics.drop([spe.replace(" ", "_")], axis=1))
    ols_df_metagenomics_X_dummies = sm.add_constant(
        ols_df_metagenomics_X_dummies)
    # OLS
    model = sm.OLS(ols_df_metagenomics_Y, ols_df_metagenomics_X_dummies).fit()
    df_tmp = pd.DataFrame({"coefs": model.params, "pvals": model.pvalues}).drop(
        'const').reset_index()
    df_tmp["species"] = [spe]*len(df_tmp)
    df_metagenomics_ols = pd.concat([df_metagenomics_ols, df_tmp])

    # Total RNA-Seq
    # Prepare df
    ols_df_totalrnaseq = df_totalrnaseq_ssu[[spe, 'trimmingscore', 'assemblytool', 'rrnasortingtool']].rename(
        columns={spe: spe.replace(" ", "_")})
    ols_df_totalrnaseq_Y = ols_df_totalrnaseq[spe.replace(" ", "_")]
    # Set up Y. Y data is not normal distributed, so we do a quantile transformation to transform into a normal distribution
    ols_df_totalrnaseq_Y = ols_df_totalrnaseq_Y.to_numpy().reshape(-1, 1)
    ols_df_totalrnaseq_Y = quantile.fit_transform(ols_df_totalrnaseq_Y)
    ols_df_totalrnaseq_Y = ols_df_totalrnaseq_Y.reshape(
        len(ols_df_totalrnaseq_Y))
    # Set up X as dummies
    ols_df_totalrnaseq_X_dummies = pd.get_dummies(
        ols_df_totalrnaseq.drop([spe.replace(" ", "_")], axis=1))
    ols_df_totalrnaseq_X_dummies = sm.add_constant(
        ols_df_totalrnaseq_X_dummies)
    # OLS
    model = sm.OLS(ols_df_totalrnaseq_Y, ols_df_totalrnaseq_X_dummies).fit()
    df_tmp = pd.DataFrame({"coefs": model.params, "pvals": model.pvalues}).drop(
        'const').reset_index()
    df_tmp["species"] = [spe]*len(df_tmp)
    df_totalrnaseq_ols = pd.concat([df_totalrnaseq_ols, df_tmp])

df_metagenomics_ols.loc[df_metagenomics_ols["pvals"]
                        <= 0.001, 'significance_cat'] = 1
df_metagenomics_ols.loc[df_metagenomics_ols["pvals"]
                        >= 0.001, 'significance_cat'] = 0.6
df_metagenomics_ols.loc[df_metagenomics_ols["pvals"]
                        >= 0.01, 'significance_cat'] = 0.35
df_metagenomics_ols.loc[df_metagenomics_ols["pvals"]
                        > 0.05, 'significance_cat'] = 0.1

df_totalrnaseq_ols.loc[df_totalrnaseq_ols["pvals"]
                       <= 0.001, 'significance_cat'] = 1
df_totalrnaseq_ols.loc[df_totalrnaseq_ols["pvals"]
                       >= 0.001, 'significance_cat'] = 0.6
df_totalrnaseq_ols.loc[df_totalrnaseq_ols["pvals"]
                       >= 0.01, 'significance_cat'] = 0.35
df_totalrnaseq_ols.loc[df_totalrnaseq_ols["pvals"]
                       > 0.05, 'significance_cat'] = 0.1

# Plot
bubble_metagenomics_ols = px.scatter(df_metagenomics_ols, x="index", y="species", size="significance_cat",
                                     color="coefs", height=550, width=840, range_color=(-df_totalrnaseq_ols["coefs"].max(), df_totalrnaseq_ols["coefs"].max()),
                                     color_continuous_scale="RdBu", template='simple_white')
# Change axes orders
bubble_metagenomics_ols.update_yaxes(categoryorder='array', categoryarray=list(reversed(['Listeria monocytogenes', 'Pseudomonas aeruginosa', 'Bacillus subtilis', 'Saccharomyces cerevisiae', 'Salmonella enterica',
                                                                                         'Escherichia coli', 'Limosilactobacillus fermentum', 'Enterococcus faecalis', 'Cryptococcus neoformans', 'Staphylococcus aureus'])))
bubble_metagenomics_ols.update_xaxes(tickangle=35, categoryorder='array', categoryarray=['trimmingscore_5', 'trimmingscore_10', 'trimmingscore_15', 'trimmingscore_20',
                                                                                         'rrnasortingtool_barrnap', 'rrnasortingtool_rrnafilter', 'rrnasortingtool_sortmerna', 'rrnasortingtool_unsorted',
                                                                                         'assemblytool_idba-tran', 'assemblytool_idba-ud', 'assemblytool_megahit', 'assemblytool_metaspades', 'assemblytool_rnaspades', 'assemblytool_spades', 'assemblytool_transabyss'])
bubble_metagenomics_ols.show()
bubble_metagenomics_ols.write_image(os.path.join(
    outdir, "ols_metagenomics.svg"), height=550, width=840)
bubble_metagenomics_ols.write_image(os.path.join(
    outdir, "ols_metagenomics.png"), height=550, width=840)

bubble_totalrnaseq_ols = px.scatter(df_totalrnaseq_ols, x="index", y="species", size="significance_cat",
                                    color="coefs", height=550, width=840, range_color=(-df_totalrnaseq_ols["coefs"].max(), df_totalrnaseq_ols["coefs"].max()),
                                    color_continuous_scale="RdBu", template='simple_white')
# Change axes orders
bubble_totalrnaseq_ols.update_yaxes(categoryorder='array', categoryarray=list(reversed(['Listeria monocytogenes', 'Pseudomonas aeruginosa', 'Bacillus subtilis', 'Saccharomyces cerevisiae', 'Salmonella enterica',
                                                                                        'Escherichia coli', 'Limosilactobacillus fermentum', 'Enterococcus faecalis', 'Cryptococcus neoformans', 'Staphylococcus aureus'])))
bubble_totalrnaseq_ols.update_xaxes(tickangle=35, categoryorder='array', categoryarray=['trimmingscore_5', 'trimmingscore_10', 'trimmingscore_15', 'trimmingscore_20',
                                                                                        'rrnasortingtool_barrnap', 'rrnasortingtool_rrnafilter', 'rrnasortingtool_sortmerna', 'rrnasortingtool_unsorted',
                                                                                        'assemblytool_idba-tran', 'assemblytool_idba-ud', 'assemblytool_megahit', 'assemblytool_metaspades', 'assemblytool_rnaspades', 'assemblytool_spades', 'assemblytool_transabyss'])
bubble_totalrnaseq_ols.show()
bubble_totalrnaseq_ols.write_image(os.path.join(
    outdir, "ols_totalrnaseq.svg"), height=550, width=840)
bubble_totalrnaseq_ols.write_image(os.path.join(
    outdir, "ols_totalrnaseq.png"), height=550, width=840)


# # Correlation testing - replaced by linear regression
# df_metagenomics_cor = pd.DataFrame()
# df_totalrnaseq_cor = pd.DataFrame()
#
# ## Drop trimming score
# df_metagenomics_genome = df_metagenomics_genome[df_metagenomics_genome['trimmingscore']=='20']
# df_totalrnaseq_genome = df_totalrnaseq_genome[df_totalrnaseq_genome['trimmingscore']=='20']
# df_metagenomics_ssu = df_metagenomics_ssu[df_metagenomics_ssu['trimmingscore']=='20']
# df_totalrnaseq_ssu = df_totalrnaseq_ssu[df_totalrnaseq_ssu['trimmingscore']=='20']
#
# ## Loop over species:
# for spe in df_metagenomics_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).columns:
#     ### Set X and Y
#     X_metagenomics = df_metagenomics_ssu[["assemblytool", "rrnasortingtool"]]
#     X_totalrnaseq = df_totalrnaseq_ssu[["assemblytool", "rrnasortingtool"]]
#     Y_metagenomics = df_metagenomics_ssu[spe]
#     Y_totalrnaseq = df_totalrnaseq_ssu[spe]
#     X_metagenomics_dummies = pd.get_dummies(X_metagenomics)
#     X_totalrnaseq_dummies = pd.get_dummies(X_totalrnaseq)
#
#     ### Metagenomics
#     cor_dic_metagenomics={}
#     for col in X_metagenomics_dummies.columns:
#         Y=Y_metagenomics.to_numpy()
#         X=np.array(X_metagenomics_dummies[col])
#         cor=stats.pointbiserialr(X,Y)
#         cor_dic_metagenomics[col]=cor
#     cor_df_metagenomics = pd.DataFrame(cor_dic_metagenomics , index=["coefficient", "p-value"]).transpose().reset_index()
#     cor_df_metagenomics = cor_df_metagenomics.drop(0)
#
#     #### Add species column for visualization
#     cor_df_metagenomics["species"] = [spe]*len(cor_df_metagenomics)
#
#     #### Add significance categories visualization
#     cor_df_metagenomics.loc[cor_df_metagenomics["p-value"] <= 0.001 , 'significance_cat'] = 1
#     cor_df_metagenomics.loc[cor_df_metagenomics["p-value"] >= 0.001 , 'significance_cat'] = 0.6
#     cor_df_metagenomics.loc[cor_df_metagenomics["p-value"] >= 0.01 , 'significance_cat'] = 0.35
#     cor_df_metagenomics.loc[cor_df_metagenomics["p-value"] > 0.05, 'significance_cat'] = 0.1
#
#     df_metagenomics_cor = pd.concat([df_metagenomics_cor, cor_df_metagenomics])
#
#     ### Total RNA-Seq
#     cor_dic_totalrnaseq={}
#     for col in X_totalrnaseq_dummies.columns:
#         Y=Y_totalrnaseq.to_numpy()
#         X=np.array(X_totalrnaseq_dummies[col])
#         cor=stats.pointbiserialr(X,Y)
#         cor_dic_totalrnaseq[col]=cor
#     cor_df_totalrnaseq = pd.DataFrame(cor_dic_totalrnaseq , index=["coefficient", "p-value"]).transpose().reset_index()
#     cor_df_totalrnaseq = cor_df_totalrnaseq.drop(0)
#
#     #### Add species column for visualization
#     cor_df_totalrnaseq["species"] = [spe]*len(cor_df_totalrnaseq)
#
#     #### Add significance categories visualization
#     cor_df_totalrnaseq.loc[cor_df_totalrnaseq["p-value"] <= 0.001 , 'significance_cat'] = 1
#     cor_df_totalrnaseq.loc[cor_df_totalrnaseq["p-value"] >= 0.001 , 'significance_cat'] = 0.6
#     cor_df_totalrnaseq.loc[cor_df_totalrnaseq["p-value"] >= 0.01 , 'significance_cat'] = 0.35
#     cor_df_totalrnaseq.loc[cor_df_totalrnaseq["p-value"] > 0.05, 'significance_cat'] = 0.1
#
#     df_totalrnaseq_cor = pd.concat([df_totalrnaseq_cor, cor_df_totalrnaseq])
#
# ## Plot
# bubble_metagenomics_cor = px.scatter(df_metagenomics_cor, x="index", y="species", size="significance_cat",
#     color="coefs", height=550, width=900, range_color=(-1,1),
#     color_continuous_scale="RdBu", template='simple_white')
# ### Change axes orders
# bubble_metagenomics_cor.update_yaxes(categoryorder='array', categoryarray= list(reversed(['Listeria monocytogenes', 'Pseudomonas aeruginosa', 'Bacillus subtilis', 'Saccharomyces cerevisiae', 'Salmonella enterica',
# 'Escherichia coli', 'Limosilactobacillus fermentum', 'Enterococcus faecalis', 'Cryptococcus neoformans', 'Staphylococcus aureus'])))
# bubble_metagenomics_cor.update_xaxes(tickangle=35, categoryorder='array', categoryarray= list(reversed(['assemblytool_idba-tran', 'assemblytool_idba-ud', 'assemblytool_megahit', 'assemblytool_metaspades', 'assemblytool_rnaspades',
# 'assemblytool_spades', 'assemblytool_transabyss', 'rrnasortingtool_barrnap', 'rrnasortingtool_rrnafilter', 'rrnasortingtool_sortmerna', 'rrnasortingtool_unsorted'])))
# bubble_metagenomics_cor.show()
# bubble_metagenomics_cor.write_image(os.path.join(outdir, "correlation_metagenomics.svg"))
# bubble_metagenomics_cor.write_image(os.path.join(outdir, "correlation_metagenomics.png"))
#
# bubble_totalrnaseq_cor = px.scatter(df_totalrnaseq_cor, x="index", y="species", size="significance_cat",
#     color="coefficient", height=550, width=700, range_color=(-1,1),
#     color_continuous_scale="RdBu", template='simple_white')
# ### Change axes orders
# bubble_totalrnaseq_cor.update_yaxes(categoryorder='array', categoryarray= list(reversed(['Listeria monocytogenes', 'Pseudomonas aeruginosa', 'Bacillus subtilis', 'Saccharomyces cerevisiae', 'Salmonella enterica',
# 'Escherichia coli', 'Limosilactobacillus fermentum', 'Enterococcus faecalis', 'Cryptococcus neoformans', 'Staphylococcus aureus'])))
# bubble_totalrnaseq_cor.update_xaxes(tickangle=35, categoryorder='array', categoryarray= list(reversed(['assemblytool_idba-tran', 'assemblytool_idba-ud', 'assemblytool_megahit', 'assemblytool_metaspades', 'assemblytool_rnaspades',
# 'assemblytool_spades', 'assemblytool_transabyss', 'rrnasortingtool_barrnap', 'rrnasortingtool_rrnafilter', 'rrnasortingtool_sortmerna', 'rrnasortingtool_unsorted'])))
# bubble_totalrnaseq_cor.show()
# bubble_totalrnaseq_cor.write_image(os.path.join(outdir, "correlation_totalrnaseq.svg"))
# bubble_totalrnaseq_cor.write_image(os.path.join(outdir, "correlation_totalrnaseq.png"))


# Coverage heatmaps
# Re-sort dfs
df_metagenomics_genome = df_metagenomics_genome.sort_values(
    ["rrnasortingtool", "assemblytool"])
df_totalrnaseq_genome = df_totalrnaseq_genome.sort_values(
    ["rrnasortingtool", "assemblytool"])
df_metagenomics_ssu = df_metagenomics_ssu.sort_values(
    ["rrnasortingtool", "assemblytool"])
df_totalrnaseq_ssu = df_totalrnaseq_ssu.sort_values(
    ["assemblytool", "rrnasortingtool"])

# Set index for readibility in plot
df_metagenomics_genome.index = df_metagenomics_genome["rrnasortingtool"] + \
    "_" + df_metagenomics_genome["assemblytool"]
df_totalrnaseq_genome.index = df_totalrnaseq_genome["rrnasortingtool"] + \
    "_" + df_totalrnaseq_genome["assemblytool"]
df_metagenomics_ssu.index = df_metagenomics_ssu["rrnasortingtool"] + \
    "_" + df_metagenomics_ssu["assemblytool"]
df_totalrnaseq_ssu.index = df_totalrnaseq_ssu["assemblytool"] + \
    "_" + df_totalrnaseq_ssu["rrnasortingtool"]

# Heatmap generation
heatmap_metagenomics_genome = px.imshow(df_metagenomics_genome.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).transpose(
), aspect="auto", width=1150, labels=dict(x="Pipeline", y="Mock community species", color="Genome coverage [%]"), zmin=0, zmax=100)
heatmap_totalrnaseq_genome = px.imshow(df_totalrnaseq_genome.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).transpose(
), aspect="auto", width=1150, labels=dict(x="Pipeline", y="Mock community species", color="Genome coverage [%]"), zmin=0, zmax=100)
heatmap_metagenomics_ssu = px.imshow(df_metagenomics_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).transpose(
), aspect="auto", width=1150, labels=dict(x="Pipeline", y="Mock community species", color="SSU coverage [%]"), zmin=0, zmax=100)
heatmap_totalrnaseq_ssu = px.imshow(df_totalrnaseq_ssu.drop(["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).transpose(
), aspect="auto", width=1150, labels=dict(x="Pipeline", y="Mock community species", color="SSU coverage [%]"), zmin=0, zmax=100)

# Plot
heatmap_metagenomics_genome.show()
heatmap_totalrnaseq_genome.show()
heatmap_metagenomics_ssu.show()
heatmap_totalrnaseq_ssu.show()

heatmap_metagenomics_genome.write_image(
    os.path.join(outdir, "coverage_metagenomics_genome.svg"))
heatmap_metagenomics_ssu.write_image(
    os.path.join(outdir, "coverage_metagenomics_ssu.svg"))
heatmap_totalrnaseq_genome.write_image(
    os.path.join(outdir, "coverage_totalrnaseq_genome.svg"))
heatmap_totalrnaseq_ssu.write_image(
    os.path.join(outdir, "coverage_totalrnaseq_ssu.svg"))

heatmap_metagenomics_genome.write_image(
    os.path.join(outdir, "coverage_metagenomics_genome.png"))
heatmap_metagenomics_ssu.write_image(
    os.path.join(outdir, "coverage_metagenomics_ssu.png"))
heatmap_totalrnaseq_genome.write_image(
    os.path.join(outdir, "coverage_totalrnaseq_genome.png"))
heatmap_totalrnaseq_ssu.write_image(
    os.path.join(outdir, "coverage_totalrnaseq_ssu.png"))


# Correlation between genome size, abundance, and coverage
genome_sizes = [2992342, 6792330, 4045677, 12071326,
                4809318, 4875441, 1905333, 2845392, 19051922, 2730326]
abundances = [94.8, 4.2, 0.7, 0.23, 0.059,
              0.058, 0.015, 0.001, 0.00015, 0.0001]
# Metagenomics
ave_coverage_metagenomics = df_metagenomics_ssu.drop(
    ["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).mean()
# ### Pearson Cor
# genome_size_pearson_metagenomics = stats.pearsonr(ave_coverage_metagenomics, genome_sizes)
# abundance_pearson_metagenomics = stats.pearsonr(ave_coverage_metagenomics, abundances)
# print("Genome size pearson metagenomics:", genome_size_pearson_metagenomics)
# print("Abundance pearson metagenomics:", abundance_pearson_metagenomics)
### Scatterplot and OLS
# vs. genome size
fig = px.scatter(pd.DataFrame({"coverage": ave_coverage_metagenomics, "genome_sizes": genome_sizes}).reset_index(),
                 x='genome_sizes', y='coverage', trendline='ols', color="index", symbol="index",
                 symbol_sequence=[0, 1, 2, 3, 4, 13, 17, 19, 20, 22],
                 color_discrete_sequence=['#f455ad', '#002684', '#B74258', '#62ecb4',
                                          '#9c5fb4', '#600008', '#5d9344', '#ff8a4e', '#01a2ed', '#AF903F'],
                 title="Metagenomics genome size", template="plotly_white", width=400, height=400)
fig.update_yaxes(range=[0, 100])
fig.update_layout(showlegend=False)
# Note: when using the color parameter above, the graph won't generate a trendline.
# Therefore, the code has to be run twice, once with the color parameter and the
# following block of code commented out and once without the color parameter and
# using the following block of code to generate the p-value for the trendline
# Get R2 and pval
# model = px.get_trendline_results(fig)
# results = model.iloc[0]["px_fit_results"]
# p_beta = results.pvalues[1]
# r_squared = results.rsquared
# pval = 'p-value = ' + '{:.5f}'.format(p_beta)
# fig.add_annotation(text=pval, x=19900000, y=95, showarrow=False, xanchor='right' , yanchor='auto')
fig.update_traces(marker_size=9)
fig.show()
fig.write_image(os.path.join(
    outdir, "scatter_metagenomics_coverage_vs_size.png"))
fig.write_image(os.path.join(
    outdir, "scatter_metagenomics_coverage_vs_size.svg"))
#### vs. abundance
fig = px.scatter(pd.DataFrame({"coverage": ave_coverage_metagenomics, "abundance": abundances}).reset_index(),
                 x='abundance', y='coverage', trendline='ols', color="index", symbol="index",
                 symbol_sequence=[0, 1, 2, 3, 4, 13, 17, 19, 20, 22],
                 color_discrete_sequence=['#f455ad', '#002684', '#B74258', '#62ecb4',
                                          '#9c5fb4', '#600008', '#5d9344', '#ff8a4e', '#01a2ed', '#AF903F'],
                 title="Metagenomics abundance", log_x=True, template="plotly_white",
                 width=400, height=400)
fig.update_yaxes(range=[0, 100])
fig.update_layout(showlegend=False)
# Note: when using the color parameter above, the graph won't generate a trendline.
# Therefore, the code has to be run twice, once with the color parameter and the
# following block of code commented out and once without the color parameter and
# using the following block of code to generate the p-value for the trendline
# Get R2 and pval
# model = px.get_trendline_results(fig)
# results = model.iloc[0]["px_fit_results"]
# p_beta = results.pvalues[1]
# r_squared = results.rsquared
# pval = 'p-value = ' + '{:.5f}'.format(p_beta)
# fig.add_annotation(text=pval, x=0.9, y=95, showarrow=False, xanchor='right' , yanchor='auto')
fig.update_traces(marker_size=9)
fig.show()
fig.write_image(os.path.join(
    outdir, "scatter_metagenomics_coverage_vs_abun.png"))
fig.write_image(os.path.join(
    outdir, "scatter_metagenomics_coverage_vs_abun.svg"))
# vs. abundance*genome size
fig = px.scatter(pd.DataFrame({"coverage": ave_coverage_metagenomics,
                               "abundanceXsize": [a*g for a, g in zip([x/100 for x in abundances], genome_sizes)]}).reset_index(),
                 x='abundanceXsize', y='coverage', trendline='ols', color="index", symbol="index",
                 symbol_sequence=[0, 1, 2, 3, 4, 13, 17, 19, 20, 22],
                 color_discrete_sequence=['#f455ad', '#002684', '#B74258', '#62ecb4',
                                          '#9c5fb4', '#600008', '#5d9344', '#ff8a4e', '#01a2ed', '#AF903F'],
                 title="Metagenomics abundanceXsize", log_x=True, template="plotly_white",
                 width=400, height=400)
fig.update_yaxes(range=[0, 100])
fig.update_layout(showlegend=False)
# Note: when using the color parameter above, the graph won't generate a trendline.
# Therefore, the code has to be run twice, once with the color parameter and the
# following block of code commented out and once without the color parameter and
# using the following block of code to generate the p-value for the trendline
# #### Get R2 and pval
# model = px.get_trendline_results(fig)
# results = model.iloc[0]["px_fit_results"]
# p_beta = results.pvalues[1]
# r_squared = results.rsquared
# pval = 'p-value = ' + '{:.5f}'.format(p_beta)
# fig.add_annotation(text=pval, x=5, y=95, showarrow=False, xanchor='right' , yanchor='auto')
fig.update_traces(marker_size=9)
fig.show()
fig.write_image(os.path.join(
    outdir, "scatter_metagenomics_coverage_vs_abunXsize.png"))
fig.write_image(os.path.join(
    outdir, "scatter_metagenomics_coverage_vs_abunXsize.svg"))

# Total RNA-Seq
ave_coverage_totalrnaseq = df_totalrnaseq_ssu.drop(
    ["assemblytool", "rrnasortingtool", "trimmingscore"], axis=1).mean()
# ### Pearson Cor
# genome_size_pearson_totalrnaseq = stats.pearsonr(ave_coverage_totalrnaseq, genome_sizes)
# abundance_pearson_totalrnaseq = stats.pearsonr(ave_coverage_totalrnaseq, abundances)
# print("Genome size pearson totalrnaseq:", genome_size_pearson_totalrnaseq)
# print("Abundance pearson totalrnaseq:", abundance_pearson_totalrnaseq)
### Scatterplot and OLS
# vs. genome size
fig = px.scatter(pd.DataFrame({"coverage": ave_coverage_totalrnaseq, "genome_sizes": genome_sizes}).reset_index(),
                 x='genome_sizes', y='coverage', trendline='ols', color="index", symbol="index",
                 symbol_sequence=[0, 1, 2, 3, 4, 13, 17, 19, 20, 22],
                 color_discrete_sequence=['#f455ad', '#002684', '#B74258', '#62ecb4',
                                          '#9c5fb4', '#600008', '#5d9344', '#ff8a4e', '#01a2ed', '#AF903F'],
                 title="Total RNA-Seq genome size", template="plotly_white", width=400, height=400)
fig.update_yaxes(range=[0, 100])
fig.update_layout(showlegend=False)
# Note: when using the color parameter above, the graph won't generate a trendline.
# Therefore, the code has to be run twice, once with the color parameter and the
# following block of code commented out and once without the color parameter and
# using the following block of code to generate the p-value for the trendline
# Get R2 and pval
# model = px.get_trendline_results(fig)
# results = model.iloc[0]["px_fit_results"]
# p_beta = results.pvalues[1]
# r_squared = results.rsquared
# pval = 'p-value = ' + '{:.5f}'.format(p_beta)
# fig.add_annotation(text=pval, x=19900000, y=95, showarrow=False, xanchor='right' , yanchor='auto')
fig.update_traces(marker_size=9)
fig.show()
fig.write_image(os.path.join(
    outdir, "scatter_totalrnaseq_coverage_vs_size.png"))
fig.write_image(os.path.join(
    outdir, "scatter_totalrnaseq_coverage_vs_size.svg"))
#### vs. abundance
fig = px.scatter(pd.DataFrame({"coverage": ave_coverage_totalrnaseq, "abundance": abundances}).reset_index(),
                 x='abundance', y='coverage', trendline='ols', color="index", symbol="index",
                 symbol_sequence=[0, 1, 2, 3, 4, 13, 17, 19, 20, 22],
                 color_discrete_sequence=['#f455ad', '#002684', '#B74258', '#62ecb4',
                                          '#9c5fb4', '#600008', '#5d9344', '#ff8a4e', '#01a2ed', '#AF903F'],
                 title="Total RNA-Seq abundance", log_x=True, template="plotly_white",
                 width=400, height=400)
fig.update_yaxes(range=[0, 100])
fig.update_layout(showlegend=False)
# Note: when using the color parameter above, the graph won't generate a trendline.
# Therefore, the code has to be run twice, once with the color parameter and the
# following block of code commented out and once without the color parameter and
# using the following block of code to generate the p-value for the trendline
# Get R2 and pval
# model = px.get_trendline_results(fig)
# results = model.iloc[0]["px_fit_results"]
# p_beta = results.pvalues[1]
# r_squared = results.rsquared
# pval = 'p-value = ' + '{:.5f}'.format(p_beta)
# fig.add_annotation(text=pval, x=0.9, y=95, showarrow=False, xanchor='right' , yanchor='auto')
fig.update_traces(marker_size=9)
fig.show()
fig.write_image(os.path.join(
    outdir, "scatter_totalrnaseq_coverage_vs_abun.png"))
fig.write_image(os.path.join(
    outdir, "scatter_totalrnaseq_coverage_vs_abun.svg"))
# vs. abundance*genome size
fig = px.scatter(pd.DataFrame({"coverage": ave_coverage_totalrnaseq,
                               "abundanceXsize": [a*g for a, g in zip([x/100 for x in abundances], genome_sizes)]}).reset_index(),
                 x='abundanceXsize', y='coverage', trendline='ols', color="index", symbol="index",
                 symbol_sequence=[0, 1, 2, 3, 4, 13, 17, 19, 20, 22],
                 color_discrete_sequence=['#f455ad', '#002684', '#B74258', '#62ecb4',
                                          '#9c5fb4', '#600008', '#5d9344', '#ff8a4e', '#01a2ed', '#AF903F'],
                 title="Total RNA-Seq abundanceXsize", log_x=True, template="plotly_white",
                 width=400, height=400)
fig.update_yaxes(range=[0, 100])
fig.update_layout(showlegend=False)
# Note: when using the color parameter above, the graph won't generate a trendline.
# Therefore, the code has to be run twice, once with the color parameter and the
# following block of code commented out and once without the color parameter and
# using the following block of code to generate the p-value for the trendline
# #### Get R2 and pval
# model = px.get_trendline_results(fig)
# results = model.iloc[0]["px_fit_results"]
# p_beta = results.pvalues[1]
# r_squared = results.rsquared
# pval = 'p-value = ' + '{:.5f}'.format(p_beta)
# fig.add_annotation(text=pval, x=5, y=95, showarrow=False, xanchor='right' , yanchor='auto')
fig.update_traces(marker_size=9)
fig.show()
fig.write_image(os.path.join(
    outdir, "scatter_totalrnaseq_coverage_vs_abunXsize.png"))
fig.write_image(os.path.join(
    outdir, "scatter_totalrnaseq_coverage_vs_abunXsize.svg"))


# # OLS on model with coverage as Y and genome size and abundance as X - note: I don't we need that as it is redundant with the code section above
# ## Metagenomics
# df = pd.DataFrame({"coverage": ave_coverage_metagenomics, "genome_sizes": genome_sizes, "abundances": abundances})
# model = sm.OLS(df['coverage'], sm.add_constant(df[['genome_sizes', 'abundances']])).fit()
# print(model.summary())
# ## Total RNA Seq
# df = pd.DataFrame({"coverage": ave_coverage_totalrnaseq, "genome_sizes": genome_sizes, "abundances": abundances})
# model = sm.OLS(df['coverage'], sm.add_constant(df[['genome_sizes', 'abundances']])).fit()
# print(model.summary())

# # 3D scatter plot of Metagenomics - note: too fancy
# margin = 0
# df = pd.DataFrame({"coverage": ave_coverage_metagenomics, "genome_sizes": genome_sizes, "abundance": abundances})
# X = df[['genome_sizes', 'abundance']]
# y = df['coverage']
# ## OLS model
# reg = LinearRegression().fit(X, y)
# # Create a mesh grid on which we will run our model
# x_min, x_max = X.genome_sizes.min() - margin, X.genome_sizes.max() + margin
# y_min, y_max = X.abundance.min() - margin, X.abundance.max() + margin
# xrange = np.arange(x_min, x_max, 1000)
# yrange = np.arange(y_min, y_max, 10)
# xx, yy = np.meshgrid(xrange, yrange)
# # Run model
# pred = reg.predict(np.c_[xx.ravel(), yy.ravel()])
# pred = pred.reshape(xx.shape)
# # Generate the plot
# fig = px.scatter_3d(df, x='genome_sizes', y='abundance', z='coverage')
# fig.update_traces(marker=dict(size=5))
# fig.add_traces(go.Surface(x=xrange, y=yrange, z=pred, name='pred_surface'))
# fig.update_layout(
#     scene = dict(zaxis = dict(range=[0,100],),),)
# fig.show()

# Aquarium plots
# Process dfs
# improve naming
df_coverage["pipeline"] = df_coverage["pipeline"].str.replace(
    "idba_ud", "idba-ud").str.replace("idba_tran", "idba-tran").str.replace("_taxonomy.txt", "")
df_taxa["pipeline"] = df_taxa["pipeline"].str.replace(
    "idba_ud", "idba-ud").str.replace("idba_tran", "idba-tran")
df_coverage["sample"] = df_coverage["sample"].str.replace(
    "DNA", "Metagenomics").str.replace("RNA", "Total RNA-Seq")
df_taxa["sample"] = df_taxa["sample"].str.replace(
    "DNA", "Metagenomics").str.replace("RNA", "Total RNA-Seq")
# split up tools used in each pipeline, samples, and sequencing types
df_coverage[["phredscore", "sortingtool", "assembler"]
            ] = df_coverage["pipeline"].str.split("_", expand=True)
df_coverage[["sample", "seqtype"]
            ] = df_coverage["sample"].str.split("_", expand=True)
df_taxa[["phredscore", "sortingtool", "assembler"]
        ] = df_taxa["pipeline"].str.split("_", expand=True)
df_taxa[["sample", "seqtype"]
        ] = df_taxa["sample"].str.split("_", expand=True)

# calculate number of prokaryotes and ratio prokaryotes:eukaryotes in taxa df
df_taxa["numprok"] = df_taxa["numarc"]+df_taxa["numbac"]
df_taxa["prok:euk"] = df_taxa["numprok"]/df_taxa["numeuk"]

# splitup dfs based on samples
df_coverage_F4 = df_coverage[df_coverage["sample"] == "F4"]
df_coverage_F5 = df_coverage[df_coverage["sample"] == "F5"]
df_coverage_F6 = df_coverage[df_coverage["sample"] == "F6"]

df_taxa_F4 = df_taxa[df_taxa["sample"] == "F4"]
df_taxa_F5 = df_taxa[df_taxa["sample"] == "F5"]
df_taxa_F6 = df_taxa[df_taxa["sample"] == "F6"]

# scatterplots coverage
F4_coverage_plot = px.scatter(df_coverage_F4, color="seqtype", x="numtaxa", y="numtaxa_80", symbol="assembler", template="plotly_white",
                              symbol_sequence=[0, 1, 2, 3, 4, 17, 18], width=530, height=400, hover_name="pipeline", title="F4",
                              labels={
                                  "seq_type": "Sequencing type",
                                  "numtaxa": "Total number of taxa",
                                  "numtaxa_80": "Number of taxa with >80% coverage",
                                  "assembler": "Assembler"
                              })
F4_coverage_plot.update_traces(marker=dict(size=8,
                                           line=dict(width=1,
                                                     color='Black')),
                               selector=dict(mode='markers'))
F4_coverage_plot.show()
F4_coverage_plot.write_image(os.path.join(
    outdir, "taxanums_F4.png"))
F4_coverage_plot.write_image(os.path.join(
    outdir, "taxanums_F4.svg"))


F5_coverage_plot = px.scatter(df_coverage_F5, color="seqtype", x="numtaxa", y="numtaxa_80", symbol="assembler", template="plotly_white",
                              symbol_sequence=[0, 1, 2, 3, 4, 17, 18], width=530, height=400, hover_name="pipeline", title="F5",
                              labels={
                                  "seq_type": "Sequencing type",
                                  "numtaxa": "Total number of taxa",
                                  "numtaxa_80": "Number of taxa with >80% coverage",
                                  "assembler": "Assembler"
                              })
F5_coverage_plot.update_traces(marker=dict(size=8,
                                           line=dict(width=1,
                                                     color='Black')),
                               selector=dict(mode='markers'))
F5_coverage_plot.show()
F5_coverage_plot.write_image(os.path.join(
    outdir, "taxanums_F5.png"))
F5_coverage_plot.write_image(os.path.join(
    outdir, "taxanums_F5.svg"))


F6_coverage_plot = px.scatter(df_coverage_F6, color="seqtype", x="numtaxa", y="numtaxa_80", symbol="assembler", template="plotly_white",
                              symbol_sequence=[0, 1, 2, 3, 4, 17, 18], width=530, height=400, hover_name="pipeline", title="F6",
                              labels={
                                  "seq_type": "Sequencing type",
                                  "numtaxa": "Total number of taxa",
                                  "numtaxa_80": "Number of taxa with >80% coverage",
                                  "assembler": "Assembler"
                              })

F6_coverage_plot.update_traces(marker=dict(size=8,
                                           line=dict(width=1,
                                                     color='Black')),
                               selector=dict(mode='markers'))
F6_coverage_plot.show()
F6_coverage_plot.write_image(os.path.join(
    outdir, "taxanums_F6.png"))
F6_coverage_plot.write_image(os.path.join(
    outdir, "taxanums_F6.svg"))

# boxplots taxa
ratio_boxplot = px.box(df_taxa, x="sample", y="prok:euk", color="seqtype",
                       template="plotly_white", hover_name="pipeline",
                       title="Ratio prokaryotes:eukaryotes", width=460, height=400,
                       labels={
                           "seqtype": "Sequencing type",
                           "prok:euk": "Ratio prokaryotes:eukaryotes",
                           "sample": "Sample",
                       })
ratio_boxplot.update_yaxes(nticks=8)
ratio_boxplot.show()
ratio_boxplot.write_image(os.path.join(
    outdir, "ratio_boxplot.png"))
ratio_boxplot.write_image(os.path.join(
    outdir, "ratio_boxplot.svg"))
