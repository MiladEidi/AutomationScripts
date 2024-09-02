#!/bin/bash
set -euo pipefail

# Variables setup
vcfFile=$1           # Input VCF file
outputPrefix=$2      # Output file prefix
refFASTA=$3          # Reference FASTA file
annovarFolder=$4     # ANNOVAR folder path
annovardbFolder=$5   # ANNOVAR databases folder
refBuild=$6          # Reference build (e.g., hg19 or hg38)

# Optional arguments for protocol and operation (default values used if not provided)
protocol="refGeneWithVer,ensGene,genomicSuperDups,clinvar_20240611,1000g2015aug_all,gnomad211_genome,gnomad_exome,exac03,kaviar,esp6500siv2_all,gme,hrcr1,iranome,abraom,intervar_20180118,dann,CADD16,revel,dbnsfp41c,avsnp150"
operation="g,g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f"

if [[ $# -ge 8 ]]; then
  protocol=$7
  operation=$8
fi

# Create output directory if it doesn't exist
mkdir -p "$(dirname "${outputPrefix}")"

# Check if required executables are available
command -v bcftools >/dev/null 2>&1 || { echo "BCFtools not found. Exiting."; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "Tabix not found. Exiting."; exit 1; }
command -v "${annovarFolder}/table_annovar.pl" >/dev/null 2>&1 || { echo "ANNOVAR script not found. Exiting."; exit 1; }

main() {
  ## Step 1: Split multiallelic variants into biallelic, perform left normalization,
  filteredVcf="${outputPrefix}.biallelic.leftnorm.vcf.gz"
  bcftools norm -m-both -f "${refFASTA}" -Oz -o "${filteredVcf}" "${vcfFile}"

  # Index the filtered VCF file
  tabix -p vcf "${filteredVcf}"

  ## Step 2: Annotate the VCF file with ANNOVAR
  vcfPlain="${outputPrefix}.unzipped.vcf"

  # Prepare VCF without genotypes and annotate using ANNOVAR
  bcftools view -G -Ou "${filteredVcf}" | bcftools annotate -x FORMAT/GT,^INFO/AC -o "${vcfPlain}"
  
  "${annovarFolder}/table_annovar.pl" "${vcfPlain}" "${annovardbFolder}" -buildver "${refBuild}" -out "${outputPrefix}.annovar" \
  -remove -protocol "${protocol}" -operation "${operation}" \
  -nastring . -vcfinput > "${outputPrefix}.annovar.log"

  ## Step 3: Adjust annotation types for consistency in downstream processing
  sed -i -e '/^##INFO=<ID=gnomAD/s/Type=String/Type=Float/' \
         -e '/^##INFO=<ID=ExAC_nontcga/s/Type=String/Type=Float/' \
         -e '/^##INFO=<ID=CADD/s/Type=String/Type=Float/' \
         -e '/^##INFO=<ID=REVEL/s/Type=String/Type=Float/' \
         -e 's/Eigen-/Eigen_/g' \
         -e 's/GERP++/GERPpp/g' \
         -e 's/PC-/PC_/g' \
         -e 's/M-CAP/M_CAP/g' \
         -e 's/fathmm-/fathmm_/g' "${outputPrefix}.annovar.${refBuild}_multianno.vcf"

  ## Step 4: Cleanup intermediate files
  rm -f "${vcfPlain}" \
        "${outputPrefix}.annovar.${refBuild}_multianno.vcf" \
        "${outputPrefix}.annovar.${refBuild}_multianno.txt" \
        "${outputPrefix}.annovar.avinput"
}

# Call the main function
main
