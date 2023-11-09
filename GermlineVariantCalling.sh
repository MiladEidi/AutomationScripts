#!/bin/bash

echo "#########AD Bioinformatics presents#########"
echo "##########www.adbioinformatics.net##########"
echo "#################2023 Nov###################"

# Define file paths
input="/mnt/g/man/Trimmed/Patient1.sam"

# Extract the file name and path
file_path=$(dirname -- "$input")
file_name=$(basename -- "$input")
file_name_no_extension="${file_name%.*}"

# Define the full path to GATK binary and databases
gatk="/mnt/g/NGSneeds/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"
fasta_Ref="/mnt/g/NGSneeds/Genome_FASTA/hg19/hg19.fa"
known_sites="/mnt/g/NGSneeds/Variant_Known_Databases/dbSNP_146_hg19/dbsnp_138.hg19.vcf"
threads=16

#Samtools commands
samtools view -b -@"$threads" -q 10 -F 4 -o "$file_path/${file_name_no_extension}.bam"  "$input" || { echo "Samtools command failed."; exit 1; }
#samtools sort -@"$threads" -m 1000M -n -o "$file_path/${file_name_no_extension}_Sorted.bam" "$file_path/${file_name_no_extension}.bam" || { echo "Samtools sorting failed."; exit 1; }
java -jar "$gatk" SortSam -I "$file_path/${file_name_no_extension}.bam" -O "$file_path/${file_name_no_extension}_Sorted.bam" -SO coordinate || { echo "SortSam sorting failed."; exit 1; }
# GATK commands
java -jar "$gatk" AddOrReplaceReadGroups -I "$file_path/${file_name_no_extension}_Sorted.bam" -O "$file_path/${file_name_no_extension}_Sorted_RG.bam" -PL Illumina -LB TwistWES -PU Unit1 -SM patient || { echo "GATK AddOrReplaceReadGroups failed."; exit 1; }
java -jar "$gatk" MarkDuplicates -I "$file_path/${file_name_no_extension}_Sorted_RG.bam" -O "$file_path/${file_name_no_extension}_Sorted_RG_MD.bam" -M "$file_path/${file_name_no_extension}_Sorted_RG_MD_Stats.txt" || { echo "GATK MarkDuplicates failed."; exit 1; }
java -jar "$gatk" BaseRecalibrator -I "$file_path/${file_name_no_extension}_Sorted_RG_MD.bam" -R "$fasta_Ref" -O "$file_path/${file_name_no_extension}_BQSRTable.txt" --known-sites "$known_sites" || { echo "GATK BaseRecalibrator failed."; exit 1; }
java -jar "$gatk" ApplyBQSR -I "$file_path/${file_name_no_extension}_Sorted_RG_MD.bam" -O "$file_path/${file_name_no_extension}_Sorted_RG_MD_BQSR.bam" -bqsr "$file_path/${file_name_no_extension}_BQSRTable.txt" || { echo "GATK ApplyBQSR failed."; exit 1; }
java -jar "$gatk" HaplotypeCaller -I "$file_path/${file_name_no_extension}_Sorted_RG_MD_BQSR.bam" -R "$fasta_Ref" -O "$file_path/${file_name_no_extension}.vcf" || { echo "GATK HaplotypeCaller failed."; exit 1; }

# Print a success message at the end
echo "Script execution completed successfully. By: Milad Eidi"
exit 0