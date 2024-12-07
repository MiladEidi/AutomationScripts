#!/bin/bash
set -euo pipefail

# Function to display usage information
usage() {
  echo "Usage: $0 <vcfFile> <outputPrefix> <refFASTA> <bed> <annovarFolder> <annovardbFolder> <refBuild> <thread> [protocol] [operation]"
  echo "  vcfFile:         Input VCF file"
  echo "  outputPrefix:    Prefix for output files"
  echo "  refFASTA:        Reference FASTA file"
  echo "  bed:             BED file for filtering off-target variants"
  echo "  annovarFolder:   Path to ANNOVAR scripts"
  echo "  annovardbFolder: Path to ANNOVAR databases"
  echo "  refBuild:        Reference build (e.g., hg19 or hg38)"
  echo "  Threads:         ANNOVAR threads"
  echo "  protocol:        (Optional) ANNOVAR protocol string"
  echo "  operation:       (Optional) ANNOVAR operation string"
  exit 1
}

# Check if enough arguments are provided
if [[ $# -lt 8 ]]; then
  echo "Error: Missing arguments."
  usage
fi

# Variables setup
vcfFile=$1           # Input VCF file
outputPrefix=$2      # Output file prefix
refFASTA=$3          # Reference FASTA file
bed=$4				 # BED file for hard filtering off-target variants
annovarFolder=$5     # ANNOVAR folder path
annovardbFolder=$6   # ANNOVAR databases folder
refBuild=$7          # Reference build (e.g., hg19 or hg38)
thread=$8

# Optional arguments for protocol and operation (default values used if not provided)
# if hg19
protocol="refGeneWithVer,knownGene,ensGene,genomicSuperDups,clinvar_20240611,1000g2015aug_all,gnomad211_genome,gnomad211_exome,exac03,kaviar,esp6500siv2_all,gme,hrcr1,iranome,abraom,intervar_20180118,dann,CADD16,revel,dbnsfp41c,avsnp150"
operation="g,g,g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f"

# if hg38
# protocol="refGeneWithVer,knownGene,ensGene,genomicSuperDups,clinvar_20240611,1000g2015aug_all,gnomad41_exome,gnomad41_genome,kaviar,esp6500siv2_all,gme,hrcr1,iranome,abraom,intervar_20180118,dann,CADD16,revel,dbnsfp41c,avsnp150"
# operation="g,g,g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f"


if [[ $# -ge 10 ]]; then
  protocol=$9
  operation=$10
fi

# Create output directory if it doesn't exist
mkdir -p "$(dirname "${outputPrefix}")"

# Check if required executables are available
command -v bcftools >/dev/null 2>&1 || { echo "BCFtools not found. Exiting."; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "Tabix not found. Exiting."; exit 1; }
command -v "${annovarFolder}/table_annovar.pl" >/dev/null 2>&1 || { echo "ANNOVAR script not found. Exiting."; exit 1; }

main() {
  ## Split multiallelic variants into biallelic, perform left normalization, and doing hard filter for off-target variants
  filteredVcf="${outputPrefix}.biallelic.leftnorm0.vcf.gz"
  bcftools norm -Oz -o "${filteredVcf}" -m-both -f "${refFASTA}" "${vcfFile}"
  
  tabix -p vcf -f "${filteredVcf}"
  
  bcftools filter -R "${bed}" -Oz -o "${outputPrefix}.biallelic.leftnorm.vcf.gz" "${filteredVcf}"

  # Index the filtered VCF file
  filteredVcf="${outputPrefix}.biallelic.leftnorm.vcf.gz"
  tabix -p vcf "${filteredVcf}"

  ## Annotate the VCF file with ANNOVAR
  vcfPlain="${outputPrefix}.unzipped.vcf"
  # Prepare VCF without genotypes and couple of info annotations to make the annotation faster as the input file would be smaller
  bcftools view -G -Ou "${filteredVcf}" | bcftools annotate -x FORMAT/GT,^INFO/AC -o "${vcfPlain}"
  
  # Annotation using ANNOVAR
  "${annovarFolder}/table_annovar.pl" "${vcfPlain}" "${annovardbFolder}" -thread "${thread}" -buildver "${refBuild}" -out "${outputPrefix}.annovar" \
  -remove -protocol "${protocol}" -operation "${operation}" \
  -nastring . -vcfinput > "${outputPrefix}.annovar.log"
  
  bgzip "${outputPrefix}.annovar.${refBuild}_multianno.vcf"
  tabix -p vcf "${outputPrefix}.annovar.${refBuild}_multianno.vcf.gz"
  
  # Now let's add the eliminated info to the annotated VCF file again
  vcfAnnotated="${outputPrefix}.annovar.GT.vcf"
  bcftools annotate -a "${filteredVcf}" --force -c INFO -Ov -o ${vcfAnnotated} "${outputPrefix}.annovar.${refBuild}_multianno.vcf.gz"

  ## Adjust annotation types for consistency in downstream processing
  sed -i -e '/^##INFO=<ID=gnomAD/s/Type=String/Type=Float/' \
         -e '/^##INFO=<ID=ExAC_nontcga/s/Type=String/Type=Float/' \
         -e '/^##INFO=<ID=CADD/s/Type=String/Type=Float/' \
         -e '/^##INFO=<ID=REVEL/s/Type=String/Type=Float/' \
	 -e '/^##INFO=<ID=DANN_score/s/Type=String/Type=Float/' \
   	 -e '/^##INFO=<ID=GME_AF/s/Type=String/Type=Float/' \
   	 -e '/^##INFO=<ID=abraom_freq/s/Type=String/Type=Float/' \
         -e '/^##INFO=<ID=DANN_score/s/Type=String/Type=Float/' \
	 -e 's/Eigen-/Eigen_/g' \
         -e 's/GERP++/GERPpp/g' \
         -e 's/PC-/PC_/g' \
         -e 's/M-CAP/M_CAP/g' \
         -e 's/fathmm-/fathmm_/g' "${outputPrefix}.annovar.GT.vcf"

  ## Compress the final file
  bgzip "${outputPrefix}.annovar.GT.vcf"
  tabix -p vcf "${outputPrefix}.annovar.GT.vcf.gz"
  
  ## Cleanup intermediate files
  rm -f "${vcfPlain}" \
  	"${outputPrefix}.biallelic.leftnorm0.vcf.gz" \
  	"${outputPrefix}.biallelic.leftnorm.vcf.gz" \
        "${outputPrefix}.annovar.${refBuild}_multianno.vcf.gz" \
        "${outputPrefix}.annovar.${refBuild}_multianno.vcf.gz.tbi"
        "${outputPrefix}.annovar.${refBuild}_multianno.txt" \
        "${outputPrefix}.annovar.avinput"
        

#BCFtools filtering command can be added here...

}

# Call the main function
main
