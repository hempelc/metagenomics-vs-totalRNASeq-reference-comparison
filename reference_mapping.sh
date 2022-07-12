#! /bin/bash

# Usage: mapping.sh refdir contigdir outdir
# refdir must contain the two folders with zymo reference genomes and ssrRNAs
# contigdir must contain all .fa contig files of pipelines
# outdir is the output directory, is generated if it does not already exists
# refdir="/Users/christopherhempel/Google Drive/PhD UoG/Project Shea/ZymoBIOMICS_refs"
# contigdir="/Users/christopherhempel/Desktop/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_MODIFIED_FINAL_FILES"
# outdir="/Users/christopherhempel/Desktop/outdir"

refdir=${1%/}
contigdir=${2%/}
outdir=${3%/}

# Define arrays to loop over
species=( "Bacillus_subtilis" "Cryptococcus_neoformans" "Enterococcus_faecalis" \
"Escherichia_coli" "Lactobacillus_fermentum" "Listeria_monocytogenes" "Pseudomonas_aeruginosa" \
"Saccharomyces_cerevisiae" "Salmonella_enterica" "Staphylococcus_aureus" )
levels=( "Genomes" "ssuRNA_consensus" )

# Make outdir if not exists
mkdir -p "${outdir}"

# Make sure no index files are in ref dirs
rm "${refdir}"/*/*.fasta.*
# Make sure no .sam/.bam files are in contigdir
rm "${contigdir}"/*.[bs]am
# Generate out file
echo -e "combo,level,species,coverage[%]" > "${outdir}"/mapping_output.csv

for level in "${levels[@]}"; do
  for species in "${species[@]}"; do
    # Pick the respective reference for each species
    ref="${refdir}"/"${level}"/$(ls "${refdir}"/"${level}"/ | grep "${species}")
    # Index it
    bwa index -p "${ref}" "${ref}"
    for contigfile in "${contigdir}"/*.fa*; do
      # Map contigs of each .fa file to the ref
      bwa mem "${ref}" "${contigfile}" > "${contigfile}".sam
      samtools sort "${contigfile}".sam > "${contigfile}".bam
      # Calculate average coverage across genome/chromosomes/plasmids
      totalbases=$(samtools coverage "${contigfile}".bam | cut -f 3 | tail -n +2 | awk '{ total += $1 } END { print total }')
      coveredbases=$(samtools coverage "${contigfile}".bam | cut -f 5 | tail -n +2 | awk '{ total += $1 } END { print total }')
      coverage=$(awk "BEGIN {print $coveredbases / $totalbases * 100}")
      echo "${contigfile##*/}","${level}","${species}","${coverage}" >> "${outdir}"/mapping_output.csv
      # Clean up generated .sam file
      rm "${contigfile}"*.[sb]am
    done
    # Clean up generated indices
    rm "${ref}".*
  done
done
