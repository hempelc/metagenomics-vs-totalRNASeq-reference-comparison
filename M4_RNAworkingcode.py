import pandas as pd
import plotly.express as px
df = pd.read_csv(r"C:\Users\sheac\ChrisProject1\M4_RNA\M4RNAmatches.csv")
df["pipeline"] = df["pipeline"].str.replace("trimmed_at_phred_", "").str.replace("_pipeline_final.txt.fasta.sam", "").str.replace("IDBA_", "IDBA-").str.replace("FIRST_", "FIRST-").str.replace("BLAST_", "BLAST-").str.lower()
df["trimming-score"] = df["pipeline"].str.split("_", n=-1, expand=True)[0]
df["rrnasorting-tool"] = df["pipeline"].str.split("_", n=-1, expand=True)[1]
df["assembly-tool"] = df["pipeline"].str.split("_", n=-1, expand=True)[2]
df["mapping-tool"] = df["pipeline"].str.split("_", n=-1, expand=True)[3]
df['database'] = df["pipeline"].str.split("_", n=-1, expand=True)[4]
df['classification-tool'] = df["pipeline"].str.split("_", n=-1, expand=True)[5]
df = df.drop(["database"], axis=1)
df = df.set_index("pipeline")
df.index.name = None
#generating proportions for each species
df["Bacillus subtilis"] = df["Bacillus subtilis"]/1558
df["Cryptococcus neoformans"] = df["Cryptococcus neoformans"]/1513
df["Enterococcus faecalis"] = df["Enterococcus faecalis"]/1562
df["Escherichia coli"] = df["Escherichia coli"]/1542
df["Lactobacillus fermentum"] = df["Lactobacillus fermentum"]/1573
df["Listeria monocytogenes"] = df["Listeria monocytogenes"]/1552
df["Pseudomonas aeruginosa"] = df["Pseudomonas aeruginosa"]/1526
df["Saccharomyces cerevisiae"] = df["Saccharomyces cerevisiae"]/1553
df["Salmonella enterica"] = df["Salmonella enterica"]/1534
df["Staphylococcus aureus"] = df["Staphylococcus aureus"]/1556
#sorting values by assembly tool
df.sort_values("assembly-tool")
#sorts rows by length of the match
df.sort_values("Bacillus subtilis")
#chris's fix to the issue I was having when first trying to generate the heatmap
pipenv install ipykernel
pip install --upgrade nbformat
#saving edited csv file above into new csv file to send to Chris for his help
df.to_csv("~/ChrisProject1/df.csv")
df
#heatmapexampleChrissent- first unsorted, then sorted by assembly tool
import pandas as pd
import plotly.express as px

df = pd.read_csv("~/ChrisProject1/df.csv")
df
#sorting pipelines by column of tools, so column A=trimming-score, B=rrnasorting-tool, C=assembly-tool, D=mapping-tool, and E=classification-tool
df = df.sort_values("trimming-score")
df = df.sort_values("rrnasorting-tool")
df = df.sort_values("assembly-tool")
df = df.sort_values("mapping-tool")
df = df.sort_values("classification-tool")

#renaming columns to include their abundances within the column name
df.rename(columns = {'Bacillus subtilis':'Bacillus subtilis-0.89%', 'Cryptococcus neoformans':'Cryptococcus neoformans-0.00089%', 'Enterococcus faecalis':'Enterococcus faecalis-0.00089%', 'Escherichia coli':'Escherichia coli-0.089%', 'Lactobacillus fermentum':'Lactobacillus fermentum-0.0089%', 'Listeria monocytogenes':'Listeria monocytogenes-89.1%', 'Pseudomonas aeruginosa':'Pseudomonas aeruginosa-8.9%', 'Saccharomyces cerevisiae':'Saccharomyces cerevisiae-0.89%', 'Salmonella enterica':'Salmonella enterica-0.089%', 'Staphylococcus aureus':'Staphylococcus aureus-0.000089%'}, inplace=True)

#sorting columns from most abundant species to least abundant species here
df = df.reindex(columns=['Listeria monocytogenes-89.1%','Pseudomonas aeruginosa-8.9%','Bacillus subtilis-0.89%','Saccharomyces cerevisiae-0.89%','Escherichia coli-0.089%','Salmonella enterica-0.089%','Lactobacillus fermentum-0.0089%','Enterococcus faecalis-0.00089%','Cryptococcus neoformans-0.00089%','Staphylococcus aureus-0.000089%'])
#this was giving me issues, removed this line of code from my own tweak to the code that chris sent me, will have to ask him what this line of code does?
df = df.iloc[: , 1:9]
#removing long line of pipeline names
df = df.reset_index().drop("index", axis=1)
#heatmap generation
fig = px.imshow(df, color_continuous_scale='peach', labels=dict(x="Species (% Abundance)", y="Pipeline number", color="Proportion"))
fig.show()
