#!/bin/bash

echo "######### Whole Exome Sequencing Variant Calling Pipeline #########"

# Change the addresses from here...
fq1_pair="/media/milad/9117284696/Course/Case1/Trimmed/S312_1P"
fq2_pair="/media/milad/9117284696/Course/Case1/Trimmed/S312_2P"

fq1_unpair="/media/milad/9117284696/Course/Case1/Trimmed/S312_1U"
fq2_unpair="/media/milad/9117284696/Course/Case1/Trimmed/S312_2U"

hisat2_dir="/usr/bin/"
hisat2_Index="/media/milad/9117284696_AD/NGSneeds/Genome_Indexes/hg19-hisat2/genome" 
gatk="/media/milad/9117284696_AD/NGSneeds/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"
fasta_Ref="/media/milad/9117284696_AD/NGSneeds/Genome_FASTA/hg19/hg19.fa"

#BQSR database
known_sites="/media/milad/9117284696_AD/NGSneeds/Variant_Known_Databases/dbSNP_146_hg19/dbsnp_138.hg19.vcf"

threads=12
bed_for_hardfiltering="/media/milad/9117284696_AD1/NGSneeds/hg19_Twist_ILMN_Exome.bed"
tmp_dir="/media/milad/9117284696_AD/"

###Remind that you also need to install BCFtools and freebayes, and tabix using "sudo apt install command"

# To here.



file_path=$(dirname -- "$fq1_pair")

echo "###############################"
echo "Step1: Aligning FASTQ reads using Hisat2 and filtering unmapped and multimapped reads while converting SAM to BAM..."
echo "###############################"
${hisat2_dir}/hisat2 -t --add-chrname -p ${threads} -x ${hisat2_Index} -1 ${fq1_pair} -2 ${fq2_pair} -U ${fq1_unpair},${fq2_unpair} | samtools view -@ ${threads} -b -q 10 -F 4 -o "$file_path/align.bam" || { echo "Hisat2 alignment failed."; exit 1; }

echo "###############################"
echo "Step2: Sorting the BAM file..."
echo "###############################"
java -jar "$gatk" SortSam --TMP_DIR $tmp_dir --VERBOSITY ERROR -I "$file_path/align.bam" -O "$file_path/align_Sorted.bam" -SO coordinate || { echo "SortSam sorting failed."; exit 1; }

echo "###############################"
echo "Step3: Adding read groups to the BAM file..."
echo "###############################"
java -jar "$gatk" AddOrReplaceReadGroups --VERBOSITY ERROR -I "$file_path/align_Sorted.bam" -O "$file_path/align_Sorted_RG.bam" -PL Illumina -LB Twist_WES -PU Unit1 -SM patient || { echo "GATK AddOrReplaceReadGroups failed."; exit 1; }

echo "###############################"
echo "Step4: Marking duplicates in the BAM file..."
echo "###############################"
java -jar "$gatk" MarkDuplicates --VERBOSITY ERROR -I "$file_path/align_Sorted_RG.bam" -O "$file_path/align_Sorted_RG_MD.bam" -M "$file_path/align_MD_Stats.txt" --REMOVE_DUPLICATES true || { echo "GATK MarkDuplicates failed."; exit 1; }

echo "###############################"
echo "Step5: Base Quality Score Recalibration (BQSR) of the BAM file..."
echo "###############################"
java -jar "$gatk" BaseRecalibrator --verbosity ERROR -I "$file_path/align_Sorted_RG_MD.bam" -R "$fasta_Ref" -O "$file_path/align_BQSRTable.txt" --known-sites "$known_sites" || { echo "GATK BaseRecalibrator failed."; exit 1; }
java -jar "$gatk" ApplyBQSR --verbosity ERROR -I "$file_path/align_Sorted_RG_MD.bam" -O "$file_path/align_Sorted_RG_MD_BQSR.bam" -bqsr "$file_path/align_BQSRTable.txt" || { echo "GATK ApplyBQSR failed."; exit 1; }

echo "###############################"
echo "Step6: Variant Calling using HaplotypeCaller..."
echo "###############################"
java -jar "$gatk" HaplotypeCaller --verbosity ERROR -I "$file_path/align_Sorted_RG_MD_BQSR.bam" -R "$fasta_Ref" -O "$file_path/Patient_GATK.vcf" || { echo "GATK HaplotypeCaller failed."; exit 1; }

echo "###############################"
echo "Step7: Variant Calling using BCFtools..."
echo "###############################"
bcftools mpileup -Ou -f "$fasta_Ref" "$file_path/align_Sorted_RG_MD_BQSR.bam" | bcftools call -mv -Ov -o "$file_path/Patient_bcftools.vcf" || { echo "BCFtools variant calling failed."; exit 1; }

echo "###############################"
echo "Step8: Variant Calling using FreeBayes..."
echo "###############################"
freebayes -f "$fasta_Ref" "$file_path/align_Sorted_RG_MD_BQSR.bam" > "$file_path/Patient_freebayes.vcf" || { echo "FreeBayes variant calling failed."; exit 1; }

echo "###############################"
echo "Step9: Merging VCF files using BCFtools..."
echo "###############################"

bgzip "$file_path/Patient_GATK.vcf"
tabix -p vcf "$file_path/Patient_GATK.vcf.gz"
bgzip "$file_path/Patient_bcftools.vcf"
tabix -p vcf "$file_path/Patient_bcftools.vcf.gz"
bgzip "$file_path/Patient_freebayes.vcf"
tabix -p vcf "$file_path/Patient_freebayes.vcf.gz"

bcftools merge --force-samples -o "$file_path/Patient_merged.vcf.gz" -Oz "$file_path/Patient_GATK.vcf.gz" "$file_path/Patient_bcftools.vcf.gz" "$file_path/Patient_freebayes.vcf.gz" || { echo "BCFtools merge failed."; exit 1; }

tabix -p vcf "$file_path/Patient_merged.vcf.gz"

echo "###############################"
echo "Step10: Hard filtering..."
echo "###############################"
java -jar "$gatk" VariantFiltration --verbosity ERROR \
   -R "$fasta_Ref" \
   -V "$file_path/Patient_merged.vcf.gz" \
   -O "$file_path/Patient_HardFiltered.vcf" \
   --filter-expression "QUAL < 20.0" --filter-name "LowQual" \
   --filter-expression "QD < 2.0" --filter-name "LowQD" \
   #--filter-expression "DP < 0" --filter-name "DepthFilter" \
   --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
   --filter-expression "MQRankSum < -12.5" --filter-name "LowMQRankSum" \
   --filter-expression "ReadPosRankSum < -8.0 || ReadPosRankSum > 8.0" --filter-name "ReadPosRankSumFilter" \
   --filter-expression "FS > 60.0" --filter-name "StrandBias" \
   --filter-expression "BaseQRankSum < -8.0" --filter-name "LowBaseQRankSum" \
   -L "$bed_for_hardfiltering" || { echo "GATK VariantFiltration failed."; exit 1; }

java -jar "$gatk" SelectVariants --verbosity ERROR \
    -R "$fasta_Ref" \
    -V "$file_path/Patient_HardFiltered.vcf" \
    -O "$file_path/Patient_HardFiltered_removed.vcf" \
    --exclude-filtered || { echo "GATK SelectVariants failed."; exit 1; }

echo "###############################"
echo "###############################"
echo "Script execution completed successfully! By: Milad Eidi"
echo "###############################"
echo "###############################"
exit 0
