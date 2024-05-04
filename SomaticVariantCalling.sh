#!/usr/bin/bash

###Some modifications need to be done on this script, be careful to use this...
###Calling variants on normal set should be done based on germline variant calling pipeline

echo "#########AD Bioinformatics presents#########"
echo "##########www.adbioinformatics.net##########"
echo "#################2023 Dec###################"

gatk='/media/ad/9117284696_AD/carrying/NGSguide/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar'
Varscan2='VarScan.v2.3.9.jar'

Thread=16
fasta_Ref='/media/ad/9117284696_AD/NGSneeds/fasta_Ref/hg19/hg19.fa'
known_sites='/media/ad/9117284696_AD/NGSneeds/Variant_Known_Databases/dbSNP_146_hg19/dbsnp_138.hg19.vcf'

Input_Normal='/media/ad/9117284696_AD/SomaticVariantCalling/Exome_Norm/Exome_Normal.bam'
Input_Tumor='/media/ad/9117284696_AD/SomaticVariantCalling/Exome_Tumor/Exome_Tumor.bam'

# Extract the files name and path
file_path_Normal=$(dirname -- "$Input_Normal")
file_name_Normal=$(basename -- "$Input_Normal")
file_name_no_extension_Normal="${file_name%.*}"
# Extract the files name and path
file_path_Tumor=$(dirname -- "$Input_Tumor")
file_name_Tumor=$(basename -- "$Input_Tumor")
file_name_no_extension_Tumor="${file_name_Tumor%.*}"

# BAM sorting
samtools sort -@ ${Thread} -o "$file_path_Normal/${file_name_no_extension_Normal}_Sorted.bam" "$file_path_Normal/$file_name_Normal"
samtools sort -@ ${Thread} -o "$file_path_Tumor/${file_name_no_extension_Tumor}_Sorted.bam" "$file_path_Tumor/$file_name_Tumor"

## Remove unsorted BAM files
## rm -f "$file_path_Normal/$file_name_Normal" "$file_path_Tumor/$file_name_Tumor"

# Read group adding
java -jar ${gatk} AddOrReplaceReadGroups -I "$file_path_Normal/${file_name_no_extension_Normal}_Sorted.bam" -O "$file_path_Normal/${file_name_no_extension_Normal}_Sorted_RG.bam" -LB WES-Somatic -PL Illumina -PU Unit1 -SM Normal
java -jar ${gatk} AddOrReplaceReadGroups -I "$file_path_Tumor/${file_name_no_extension_Tumor}_Sorted.bam" -O "$file_path_Tumor/${file_name_no_extension_Tumor}_Sorted_RG.bam" -LB WES-Somatic -PL Illumina -PU Unit1 -SM Tumor

## Remove intermediate files
## rm -f "$file_path_Normal/${file_name_no_extension_Normal}_Sorted.bam" "$file_path_Tumor/${file_name_no_extension_Tumor}_Sorted.bam"

# Mark duplicated reads
java -jar ${gatk} MarkDuplicates -I "$file_path_Normal/${file_name_no_extension_Normal}_Sorted_RG.bam" -O "$file_path_Normal/${file_name_no_extension_Normal}_Sorted_RG_MD.bam" -M "$file_path_Normal/${file_name_no_extension_Normal}.stat"
java -jar ${gatk} MarkDuplicates -I "$file_path_Tumor/${file_name_no_extension_Tumor}_Sorted_RG.bam" -O "$file_path_Tumor/${file_name_no_extension_Tumor}_Sorted_RG_MD.bam" -M "$file_path_Tumor/${file_name_no_extension_Tumor}.stat"

## Remove intermediate files
## rm -f "$file_path_Normal/${file_name_no_extension_Normal}_Sorted_RG.bam" "$file_path_Tumor/${file_name_no_extension_Tumor}_Sorted_RG.bam"

# Base Quality Score Recalibration
for input in "$file_path_Normal/${file_name_no_extension_Normal}_Sorted_RG_MD.bam" "$file_path_Tumor/${file_name_no_extension_Tumor}_Sorted_RG_MD.bam"; do
    java -jar ${gatk} BaseRecalibrator -I ${input} -O ${input}_bqsr.table -R ${fasta_Ref} --known-sites ${known_sites}
    java -jar ${gatk} ApplyBQSR -bqsr ${input}_bqsr.table -I ${input} -O ${input}_BQSR.bam
done

# Varscan2 somatic variant caller
samtools mpileup --no-BAQ -f ${fasta_Ref} ${Input_Normal}_Sorted_RG_MD_BQSR.bam ${Input_Tumor}_Sorted_RG_MD_BQSR.bam | java -jar ${Varscan2} somatic - "$file_path_Normal/Varscan2.vcf" --mpileup 1 --output-vcf

## Should be run by your own
# java -jar VarScan.v2.3.9.jar processSomatic '/media/ad/9117284696_AD/SomaticVariantCalling/Varscan/exome.snp.vcf'
# java -jar VarScan.v2.3.9.jar processSomatic '/media/ad/9117284696_AD/SomaticVariantCalling/Varscan/exome.indel.vcf'


echo "
Successfully Finished!
Presented to you by
www.ADBioinformatics.net
2023 Dec"
