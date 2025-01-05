#!/bin/bash

# Define variables
sample="/media/milad/9117284696_AD/SamaneFASTQ/328/HSR000328.120912.lane3.trimmed"

gatk="/media/milad/9117284696_AD/NGSneeds/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"
fasta="/media/milad/9117284696_AD/NGSneeds/Genome_Indexes/hisat_grch38/GRCh38.p14.genome.fa"
knownsites="/media/milad/9117284696_AD/NGSneeds/Variant_Known_Databases/dbSNP_146_hg38/dbsnp_146.hg38.vcf"
sample_name=$(basename "$sample" "_trimmed.Aligned.sortedByCoord.out")
hisat2_Index="/media/milad/9117284696_AD/NGSneeds/Genome_Indexes/hisat_grch38/genome"
tmp_dir="/media/milad/9117284696_AD/Softwares/"

# Main Steps
date
echo "Alignment of the reads..."
hisat2 -t -p 12 -x ${hisat2_Index} -U $sample.fastq.gz | samtools view -@ 12 -b -q 10 -F 4 -o "$sample.bam" || { echo "Hisat2 alignment failed."; exit 1; }

date
echo "Sorting alignment..."
java -jar "$gatk" SortSam --TMP_DIR $tmp_dir --VERBOSITY ERROR -I "$sample.bam" -O "$sample.sort.bam" -SO coordinate || { echo "SortSam sorting failed."; exit 1; }

date
echo "Adding Read Groups..."
java -jar "$gatk" AddOrReplaceReadGroups \
    -I "$sample.sort.bam" \
    -O "$sample.rg.bam" \
    -LB RNA-seq \
    -PL illumina \
    -SM "$sample_name" \
    -PU 10034 \
    --CREATE_INDEX true || { echo "Error adding read groups"; exit 1; }

date
echo "Marking Duplicates..."
java -jar "$gatk" MarkDuplicates \
    --VERBOSITY ERROR \
    -I "$sample.rg.bam" \
    -O "$sample.RG.MD.bam" \
    -M "$sample.MD_Stats.txt" \
    --REMOVE_DUPLICATES true || { echo "Error in MarkDuplicates"; exit 1; }

rm "$sample.rg.bam" || echo "Warning: Failed to remove intermediate files."

date
echo "Running SplitNCigarReads..."
java -jar "$gatk" SplitNCigarReads \
    -R "$fasta" \
    -I "$sample.RG.MD.bam" \
    -O "$sample.RG.MD.cigar.bam" \
    --verbosity ERROR || { echo "Error in SplitNCigarReads"; exit 1; }

rm "$sample.RG.MD.bam" || echo "Warning: Failed to remove intermediate files."

date
echo "Base Quality Score Recalibration..."
java -jar "$gatk" BaseRecalibrator \
    -R "$fasta" \
    -I "$sample.RG.MD.cigar.bam" \
    --use-original-qualities \
    --known-sites "$knownsites" \
    -O "$sample.recal.table" \
    --verbosity ERROR || { echo "Error in BaseRecalibrator"; exit 1; }

date
echo "Applying BQSR..."
java -jar "$gatk" ApplyBQSR \
    --add-output-sam-program-record \
    -R "$fasta" \
    -I "$sample.RG.MD.cigar.bam" \
    --use-original-qualities \
    -O "$sample.RG.MD.cigar.recal.bam" \
    --bqsr-recal-file "$sample.recal.table" \
    --verbosity ERROR || { echo "Error applying BQSR"; exit 1; }

rm "$sample.RG.MD.cigar.bam" || echo "Warning: Failed to remove intermediate files."
samtools index "$sample.RG.MD.cigar.recal.bam"

date
echo "Calling Variants with HaplotypeCaller..."
java -jar "$gatk" HaplotypeCaller \
    -R "$fasta" \
    -I "$sample.RG.MD.cigar.recal.bam" \
    --standard-min-confidence-threshold-for-calling 20 \
    --dont-use-soft-clipped-bases \
    -O "$sample.hc.vcf" \
    --verbosity ERROR || { echo "Error in HaplotypeCaller"; exit 1; }

date
echo "Pipeline completed successfully."

