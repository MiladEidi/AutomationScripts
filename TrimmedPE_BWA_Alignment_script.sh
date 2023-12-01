#!/usr/bin/bash

echo "Paired-end alignment started..."

fq1_pair=SRR13963617_ForwardPairs.fastq.gz
fq2_pair=SRR13963617_ReversePairs.fastq.gz
fq3_unpair=SRR13963617_ForwardUnpairs.fastq.gz
fq4_unpair=SRR13963617_ReverseUnpairs.fastq.gz

SAM1=SRR13963617_BWA.sam
SAM2=SRR13963617_BWA_FUnpair.sam
SAM3=SRR13963617_BWA_RUnpair.sam

bwa_index='/media/ad/9117284696_AD1/NGSneeds/Genome_Indexes/GRCh37_BWAIndex/genome.fa'

Thread_Numbers=12


bwa mem -t ${Thread_Numbers} -o ${SAM1} ${bwa_index} ${fq1_pair} ${fq2_pair}
bwa mem -t ${Thread_Numbers} -o ${SAM2} ${bwa_index} ${fq3_unpair}
bwa mem -t ${Thread_Numbers} -o ${SAM3} ${bwa_index} ${fq4_unpair}

samtools sort -@ ${Thread_Numbers} -o ${SAM1}_Sorted.sam ${SAM1}
samtools sort -@ ${Thread_Numbers} -o ${SAM2}_Sorted.sam ${SAM2}
samtools sort -@ ${Thread_Numbers} -o ${SAM3}_Sorted.sam ${SAM3}

rm -f ${SAM1}
rm -f ${SAM2}
rm -f ${SAM3} 

samtools merge -@ ${Thread_Numbers} -O SAM ${SAM1}MergedSortedBWA.sam ${SAM1}_Sorted.sam ${SAM2}_Sorted.sam ${SAM3}_Sorted.sam

rm -f ${SAM1}_Sorted.sam
rm -f ${SAM2}_Sorted.sam
rm -f ${SAM3}_Sorted.sam



echo "



Paired-End alignment finished successfully!  
By
Milad Eidi
AD Bioinformatics
Feb 2023
"