#!/bin/bash

# Version 0.1
# Written by Natalie Wright (nwrigh06@uoguelph.ca) and Chris Hempel (hempelc@uoguelph.ca)
# 6 Jul 2022

# This is a modified pipeline for Chris Hempel's second PhD chapter, based on
# the pipelines from his first chapter.
# You can find the original pipeline here:
# https://github.com/hempelc/metagenomics-vs-totalRNASeq
# Only the Trimming, rRNAFiltering, and Assembly step of the pipelines are
# performed

# It trims raw paired-end input sequences at 4 PHRED scores, filters rRNA
# with 4 approaches, and uses 7 assemblers.

# The output is a folder called METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_MODIFIED/FINAL_FILES/
# that contains .fa files for all combinations of the above described steps.

# The pipeline requires the following subscripts, which are all located in the
# subscripts/ directory:
	# fasta_to_tab, fastqc_on_R1_R2_and_optional_trimming.sh, deinterleave_fastq_reads.sh,

# IMPORTANT: you have to manually specify the location of Trimmomatic
trimmomatic="/hdd1/programs_for_pilot/Trimmomatic-0.39/trimmomatic-0.39.jar"

# The pipeline requires the following programs/python packages (versions we used
# when writing this script are indicated in brackets):
	# FastQC (0.11.5), Trimmomatic (0.39), sortmeRNA (4.0.0), barrnap (0.9),
	# rRNAFILTER (1.1), SPADES (3.14.0)[note: runs with the --meta and --rna
	# options for METASPADES and RNASPADES], MEGAHIT (1.2.9), IDBA-UD (1.1.1),
	# IDBA_tran (1.1.1), Trinity (2.10.0),	Trans-ABySS (2.0.1), seqtk (1.2-r94),
	# python module ete3 (3.1.2)

	# Note 1: we had to edit IDBA prior to compiling it because it didn't work
	# using long reads and the -l option. This seems to be a common problem and
	# can be circumvented following for example the instructions in
	# http://seqanswers.com/forums/showthread.php?t=29109, and see also
	# https://github.com/loneknightpy/idba/issues/26

	# Note 2: we had to install ete3 via conda in a conda environment, so we have
	# to activate that environment when running this script:
	# Based on https://github.com/conda/conda/issues/7980
	source /hdd1/chempel/programs/miniconda/etc/profile.d/conda.sh
	conda activate ete3 # ete3 is our conda environment in which we installed ete3


cmd="$0 $@" # Make variable containing the entire entered command to print command to logfile
usage="$(basename "$0") -1 <R1.fastq> -2 <R2.fastq> [-t <n>]
Usage:
	-1 Full path to forward reads in .fastq/.fq format
	-2 Full path to reverse reads in .fastq/.fq format
	-t Number of threads (default:16)
	-h Display this help and exit"

# Set default options:
threads='16'

# Set specified options:
while getopts ':1:2:t:h' opt; do
 	case "${opt}" in
		1) forward_reads="${OPTARG}" ;;
		2) reverse_reads="${OPTARG}" ;;
		t) threads="${OPTARG}" ;;
		h) echo "$usage"
			 exit ;;
		:) printf "Option -$OPTARG requires an argument."
			 echo -e "\n$usage"
			 exit ;;
		\?) printf "Invalid option: -$OPTARG"
		   echo -e "\n$usage"
		   exit
	  esac
done
shift $((OPTIND - 1))

# Check if required options are set:
if [[ -z "$forward_reads" || -z "$reverse_reads" ]]
then
   echo -e "-1 and -2 must be set.\n"
   echo -e "$usage\n\n"
   echo -e "Exiting script.\n"
   exit
fi

# Define functions to print steps with time
start=$(date +%s)

step_description_and_time_first () {
	echo -e "\n++++++++ [$(date +%H:%M:%S)] ${1} [Runtime: $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ++++++++\n"
}
step_description_and_time_second () {
	echo -e "\n======== [$(date +%H:%M:%S)] ${1} [Runtime: $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m] ========\n" #" adding outcommented quote here to fix bug in colouring scheme of personal text editor
}

##################### Write start time and options to output ######################


# Make open bracket to later tell script to write everything that follows into a logfile:
(

# Define starting time of script for total runtime calculation:
start=$(date +%s)
echo -e "\n\nSTART RUNNING SCRIPT AT $(date)\n"
echo -e "=================================================================\n\n"


# Output specified options:
step_description_and_time_second "OPTIONS"

echo -e "Forward reads were defined as $forward_reads.\n"
echo -e "Reverse reads were defined as $reverse_reads.\n"
echo -e "Number of threads was set to $threads.\n"
echo -e "Script started with full command: $cmd\n"



######################### Start of the actual script ################################
step_description_and_time_first "START RUNNING SCRIPT"


# Make output directory and directory for final files:
mkdir METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_MODIFIED/
mkdir METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_MODIFIED/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_MODIFIED_FINAL_FILES/
cd METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_MODIFIED/

# Save full current path in variable to make navigation between directories easier:
base_directory=$(pwd)

######################### Step 1: trimming ################################
step_description_and_time_first "START STEP 1: TRIMMING AND ERROR CORRECTION"

# Trimming is done with a separate subscript:
fastqc_on_R1_R2_and_optional_trimming.sh \
-T $trimmomatic \
-1 $forward_reads -2 $reverse_reads -t yes -p $threads
mv fastqc_on_R1_R2_and_optional_trimming_output step_1_trimming/

# # Running error correction module of SPAdes on all trimmed reads
for trimming_results in step_1_trimming/trimmomatic/*; do
	trim_phred=$(echo ${trimming_results##*/} | sed 's/trimmed_at_phred_//g' | sed 's/_.*$//g')
	step_description_and_time_second "ERROR-CORRECTING READS IN FOLDER $trimming_results"
	spades.py -1 $trimming_results/*1P.fastq -2 $trimming_results/*2P.fastq \
	--only-error-correction --disable-gzip-output -o $trimming_results/error_correction \
	-t $threads
	mv $trimming_results/error_correction/corrected/*1P*.fastq \
	$trimming_results/error_correction/corrected/*2P*.fastq $trimming_results
	# Rename weird name of error-corrected reads:
	R1=$(echo $trimming_results/*1P.00.0_0.cor.fastq) \
  && 	mv $trimming_results/*1P.00.0_0.cor.fastq ${R1%.00.0_0.cor.fastq}_error_corrected.fastq
  sed -r -i 's/ BH:.{2,6}//g' ${R1%.00.0_0.cor.fastq}_error_corrected.fastq
	R2=$(echo $trimming_results/*2P.00.0_0.cor.fastq) \
  && mv $trimming_results/*2P.00.0_0.cor.fastq ${R2%.00.0_0.cor.fastq}_error_corrected.fastq
  sed -r -i 's/ BH:.{2,6}//g' ${R2%.00.0_0.cor.fastq}_error_corrected.fastq
	rm -r $trimming_results/error_correction/
	step_description_and_time_second "FINISHED ERROR-CORRECTING READS IN FOLDER $trimming_results"
done

step_description_and_time_first "FINISHED STEP 1: TRIMMING AND ERROR CORRECTION"


######################### Step 2: rRNA sorting ################################

# NOTE: To enable each combination of programs in each step, we run nested for
# loops and close them all at the very end of the script.

# For loop 1: loop over the folders for trimmed data at different PHRED scores
# for rRNA filtration:
for trimming_results in step_1_trimming/trimmomatic/*; do
	mkdir $trimming_results/step_2_rrna_sorting/
	cd $trimming_results/step_2_rrna_sorting/
	trim_phred=$(echo ${trimming_results##*/} | sed 's/trimmed_at_phred_//g' | sed 's/_.*$//g')

	step_description_and_time_first "START STEP 2: rRNA SORTING OF TRIMMED READS IN FOLDER $trimming_results"

	step_description_and_time_second "CONVERT READS IN FASTA FORMAT FOR rRNAFILTER AND BARRNAP"
	mkdir reads_in_fasta_format/
	fq2fa ../*1P_error_corrected.fastq reads_in_fasta_format/R1.fa
	fq2fa ../*2P_error_corrected.fastq reads_in_fasta_format/R2.fa
	step_description_and_time_second "READS TO FASTA CONVERSION DONE"

	step_description_and_time_second "RUNNING SORTMERNA"
	mkdir SORTMERNA/
	sortmerna --ref /hdd2/databases/sortmerna_silva_databases/silva-bac-16s-id90.fasta \
	--ref /hdd2/databases/sortmerna_silva_databases/silva-arc-16s-id95.fasta \
	--ref /hdd2/databases/sortmerna_silva_databases/silva-euk-18s-id95.fasta \
  --ref /hdd2/databases/sortmerna_silva_databases/silva-euk-28s-id98.fasta \
  --ref /hdd2/databases/sortmerna_silva_databases/silva-arc-23s-id98.fasta \
  --ref /hdd2/databases/sortmerna_silva_databases/silva-bac-23s-id98.fasta \
  --ref /hdd2/databases/sortmerna_silva_databases/rfam-5.8s-database-id98.fasta \
  --ref /hdd2/databases/sortmerna_silva_databases/rfam-5s-database-id98.fasta \
	--reads ../*1P_error_corrected.fastq --reads ../*2P_error_corrected.fastq \
	--paired_in --out2 -other -fastx 1 -num_alignments 1 -v -workdir SORTMERNA/ \
	--threads 1:1:$threads
	step_description_and_time_second "SORTMERNA DONE"

	step_description_and_time_second "RUNNING rRNAFILTER"
	mkdir rRNAFILTER/
	cd rRNAFILTER/
	# rRNAFilter only worked for us when we started it within the directory
	# containing the .jar file. To simplify switching to that directory, we simply
	# download the small program within the script and delete it after usage:
	wget http://hulab.ucf.edu/research/projects/rRNAFilter/software/rRNAFilter.zip
	unzip rRNAFilter.zip
	cd rRNAFilter/
	# We use 7GB for the rRNAFilter .jar, as shown in the rRNAFilter manual:
	java -jar -Xmx7g rRNAFilter_commandline.jar \
	-i ../../reads_in_fasta_format/R1.fa -r 0
	java -jar -Xmx7g rRNAFilter_commandline.jar \
	-i ../../reads_in_fasta_format/R2.fa -r 0
	mv ../../reads_in_fasta_format/R*.fa_rRNA ..
	cd ..
	rm -r rRNAFilter rRNAFilter.zip
	# We want to keep paired reads, so we extract all rRNA read names that were
	# found in R1 and R2, save them as one list, and extract all reads from both
	# R1 and R2 reads. That way, even if only one read from a pair was identified
	# as rRNA, we keep the pair of reads:
	fasta_to_tab R1.fa_rRNA | cut -f 1 | cut -f1 -d " " > names.txt
	fasta_to_tab R2.fa_rRNA | cut -f 1 | cut -f1 -d " " >> names.txt
	sort -u names.txt > names_sorted.txt
	seqtk subseq ../reads_in_fasta_format/R1.fa names_sorted.txt \
	> rRNAFilter_paired_R1.fa
	seqtk subseq ../reads_in_fasta_format/R2.fa names_sorted.txt \
	> rRNAFilter_paired_R2.fa
	rm names_sorted.txt names.txt
	cd ..
	step_description_and_time_second "rRNAFILTER DONE"

	step_description_and_time_second "RUNNING BARRNAP"
	mkdir BARRNAP/
	for kingdom in euk bac arc; do # barrnap needs to be run on kingdoms separately
		step_description_and_time_second "RUNNING BARRNAP ON KINGDOM $kingdom AND R1 READS"
		barrnap --quiet --lencutoff 0.000001 --reject 0.000001 --kingdom $kingdom \
		--threads $threads --outseq BARRNAP/${kingdom}_reads1.fa \
		reads_in_fasta_format/R1.fa
		step_description_and_time_second "RUNNING BARRNAP ON KINGDOM $kingdom AND R2 READS"
		barrnap --quiet --lencutoff 0.000001 --reject 0.000001 --kingdom $kingdom \
		--threads $threads --outseq BARRNAP/${kingdom}_reads2.fa \
		reads_in_fasta_format/R2.fa
		rm reads_in_fasta_format/*.fai
		sed 's/.*::/>/g' BARRNAP/${kingdom}_reads1.fa | sed 's/:[^:]*$//g' \
		> BARRNAP/${kingdom}_reads1_edited.fa
		sed 's/.*::/>/g' BARRNAP/${kingdom}_reads2.fa | sed 's/:[^:]*$//g' \
		> BARRNAP/${kingdom}_reads2_edited.fa
	done
	# Concatenating results from the three kingdoms and R1 and R2 files
	cat BARRNAP/*edited.fa > BARRNAP/all_reads.fa
	# We want to keep paired reads, so we extract all rRNA read names that were
	# found in R1 and R2 for the three kingdoms (in all_reads.fa), save them as
	# one list, and extract all reads from both R1 and R2 reads. That way, even if
	# only one read from a pair was identified as rRNA, we keep the pair of reads:
	fasta_to_tab BARRNAP/all_reads.fa | cut -f 1 | cut -f1 -d " " | sort -u \
	> BARRNAP/names_sorted.txt
	seqtk subseq reads_in_fasta_format/R1.fa BARRNAP/names_sorted.txt \
	> BARRNAP/barrnap_paired_R1.fa
	seqtk subseq reads_in_fasta_format/R2.fa BARRNAP/names_sorted.txt \
	> BARRNAP/barrnap_paired_R2.fa
	rm BARRNAP/names_sorted.txt
	step_description_and_time_second "BARRNAP DONE"

	step_description_and_time_second "MAKING FOLDER UNSORTED/ AND COPYING UNSORTED READS IN THERE TO KEEP THE FOLDER STRUCTURE CONSTANT"
	mkdir UNSORTED/
	cp ../*1P_error_corrected.fastq ../*2P_error_corrected.fastq UNSORTED/

	step_description_and_time_first "FINISHED STEP 2: rRNA SORTING OF TRIMMED READS IN FOLDER $trimming_results"


	######################### Step 3: Assembly ################################

	# For loop 2: loop over the folders for rRNA filtered data for assembly:
  for rrna_filter_results in rRNAFILTER SORTMERNA BARRNAP UNSORTED; do
		if [[ $rrna_filter_results == 'rRNAFILTER' ]]; then
			R1_sorted='rRNAFILTER/rRNAFilter_paired_R1.fa'
			R2_sorted='rRNAFILTER/rRNAFilter_paired_R2.fa'
  	elif [[ $rrna_filter_results == 'SORTMERNA' ]]; then
			R1_sorted='SORTMERNA/out/aligned_fwd.fastq'
			R2_sorted='SORTMERNA/out/aligned_rev.fastq'
		elif [[ $rrna_filter_results == 'BARRNAP' ]]; then
			R1_sorted='BARRNAP/barrnap_paired_R1.fa'
			R2_sorted='BARRNAP/barrnap_paired_R2.fa'
	  else # UNSORTED
			R1_sorted='UNSORTED/*1P_error_corrected.fastq'
			R2_sorted='UNSORTED/*2P_error_corrected.fastq'
		fi

		step_description_and_time_first "START STEP 3: ASSEMBLY OF TRIMMED READS IN FOLDER $trimming_results/ AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/"
		mkdir $rrna_filter_results/step_3_assembly/
		cd $rrna_filter_results/step_3_assembly/

		step_description_and_time_second "RUNNING SPADES"
		mkdir SPADES/
		spades.py -1 ../../$R1_sorted -2 ../../$R2_sorted --only-assembler \
		-o SPADES/ -t $threads
		step_description_and_time_second "SPADES DONE"

		step_description_and_time_second "RUNNING METASPADES"
		mkdir METASPADES/
		spades.py --meta -1 ../../$R1_sorted -2 ../../$R2_sorted --only-assembler \
		-o METASPADES/ -t $threads
		step_description_and_time_second "METASPADES DONE"

		step_description_and_time_second "RUNNING MEGAHIT"
		megahit --presets meta-large -t $threads -1 ../../$R1_sorted \
		-2 ../../$R2_sorted -o MEGAHIT/
		step_description_and_time_second "MEGAHIT DONE"

		step_description_and_time_second "RUNNING IDBA_UD"
		# Note: we had to edit IDBA prior to compiling it because it didn't work
		# using long reads and the -l option. This seems to be a common problem and
		# can be circumvented following for example the instructions in
		# http://seqanswers.com/forums/showthread.php?t=29109, and see also
		# https://github.com/loneknightpy/idba/issues/26
		# IDBA_UD only takes interleaved fasta files
		fq2fa --merge --filter ../../$R1_sorted ../../$R2_sorted idba_ud_input.fa
		idba_ud --num_threads $threads -r idba_ud_input.fa -o IDBA_UD/
		mv idba_ud_input.fa IDBA_UD/
		step_description_and_time_second "IDBA_UD DONE"

		step_description_and_time_second "RUNNING RNASPADES"
		mkdir RNASPADES/
		spades.py --rna -1 ../../$R1_sorted -2 ../../$R2_sorted --only-assembler \
		-o RNASPADES/ -t $threads
		step_description_and_time_second "RNASPADES DONE"

		step_description_and_time_second "RUNNING IDBA_TRAN"
		# IDBA_TRAN only takes interleaved fasta files
		fq2fa --merge ../../$R1_sorted ../../$R2_sorted idba_tran_input.fa
		idba_tran --num_threads $threads -l idba_tran_input.fa -o IDBA_TRAN/
		mv idba_tran_input.fa IDBA_TRAN/
		step_description_and_time_second "IDBA_TRAN DONE"

		# step_description_and_time_second "RUNNING TRINITY"
		# # Barrnap and rRNAFilter output fasta files which has to be indicated to Trinity:
    # if [[ $rrna_filter_results == "rRNAFILTER" \
		# || $rrna_filter_results == "BARRNAP" ]]; then
  	# 	Trinity --seqType fa --max_memory 20G --left ../../$R1_sorted --right \
    #   ../../$R2_sorted --CPU $threads --output TRINITY/
    # else
    #   Trinity --seqType fq --max_memory 20G --left ../../$R1_sorted --right \
    #   ../../$R2_sorted --CPU $threads --output TRINITY/
    # fi
    # cat TRINITY/Trinity.fasta | sed 's/ len/_len/g' \
		# > TRINITY/Trinity_with_length.fasta  # Edit for universal format
		# step_description_and_time_second "TRINITY DONE"

		step_description_and_time_second "RUNNING TRANSABYSS"
    transabyss --pe ../../$R1_sorted ../../$R2_sorted --threads $threads \
    --outdir TRANSABYSS/
    sed 's/ /_/g' TRANSABYSS/transabyss-final.fa \
		> TRANSABYSS/transabyss-final_edited.fa # Edit for universal format
		step_description_and_time_second "TRANSABYSS DONE"

		step_description_and_time_first "FINISHED STEP 3: ASSEMBLY OF TRIMMED READS IN FOLDER $trimming_results/ AND rRNA FILTERED READS IN FOLDER $rrna_filter_results/"


		######################### Step 4: Mapping ################################

		# For loop 3: loop over the folders for assembled data for mapping:
		for assembly_results in SPADES METASPADES MEGAHIT IDBA_UD RNASPADES IDBA_TRAN TRANSABYSS; do

			if [[ $assembly_results == 'SPADES' ]]; then
				scaffolds='SPADES/scaffolds.fasta'
			elif [[ $assembly_results == 'METASPADES' ]]; then
				scaffolds='METASPADES/scaffolds.fasta'
			elif [[ $assembly_results == 'MEGAHIT' ]]; then
				scaffolds='MEGAHIT/final.contigs.fa'
			elif [[ $assembly_results == 'IDBA_UD' ]]; then
				scaffolds='IDBA_UD/scaffold.fa'
			elif [[ $assembly_results == 'RNASPADES' ]]; then
				scaffolds='RNASPADES/transcripts.fasta'
			elif [[ $assembly_results == 'IDBA_TRAN' ]]; then
				scaffolds='IDBA_TRAN/contig.fa'
			elif [[ $assembly_results == 'TRINITY' ]]; then
				scaffolds='TRINITY/Trinity_with_length.fasta'
			else # TRANSABYSS
				scaffolds='TRANSABYSS/transabyss-final_edited.fa'
			fi

			cp $scaffolds ${base_directory}/METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_MODIFIED_FINAL_FILES/trimmed_at_phred_${trim_phred}_${rrna_filter_results}_${assembly_results}.fa

				# After each loop, we have to go back to the directory we started the
				# loop in, and we're using realtive paths for that, making use of the
				# variables generated in the previous loops
			cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/${rrna_filter_results}/step_3_assembly/)
		done
		cd $(realpath --relative-to=$(pwd) ${base_directory}/${trimming_results}/step_2_rrna_sorting/)
	done
	cd $(realpath --relative-to=$(pwd) $base_directory)
done

# Display runtime
echo -e "=================================================================\n"
echo "SCRIPT DONE AFTER $((($(date +%s)-$start)/3600))h $(((($(date +%s)-$start)%3600)/60))m"

# Write output to console and log file
) 2>&1 | tee METAGENOMICS_METATRANSCRIPTOMICS_PIPELINE_MODIFIED_LOG.txt
