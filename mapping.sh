
folder1 = folder_with_ssu_references

folder2 = folder_with_genome_references

index both:

  bwa index -p test_$file Pseudomonas_aeruginosa_16S_170923.fasta

  for file in *.fasta; do bwa index -p test_$file $file; done

Turn contigs info into fasta:
for file in *.txt; do tail -n +2 $file | awk '{print $1"\t"$NF}'| sed 's/^/>/g' | sed 's/\t/\n/g' > $file.fasta ;done

mapping into sam:
for file in *.txt.fasta; do bwa mem test_indexreferencefilename $file > $file.sam; done

grep longest match:
for i in *.txt.fasta.sam; do longest_match=$(grep -v ^@ $i | cut -f 14 | grep AS | sed 's/AS:i://g' | sort -n | tail -n 1); echo $i,$longest_match >> output.csv; done
