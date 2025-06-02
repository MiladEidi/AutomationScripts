#!/usr/bin/env Rscript

# ---------------------------- #
#  Automated Variant Filtering #
# ---------------------------- #

args <- commandArgs(trailingOnly = TRUE)

# Input validation
if (length(args) != 1) {
  stop("Usage: Rscript filter_variants.R <tab-delimitted annovar annotated variants>")
}

input_file <- args[1]
cat("✅ Filtering Started...\n📄 Input file:", input_file, "\n")
# Load Data
Variants <- read.delim(file = "/media/milad/9117284696_AD/Tohid/Milad/L6632.annotated.hg38_multianno.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("✅ Reading input file finished!\n")
# Helper to safely convert to numeric and filter based on threshold
filter_AF_less <- function(df, colname, threshold = 0.05) {
    if (!colname %in% names(df)) {
      warning(paste("Column", colname, "not found in the data frame. Function skipped."))
      return(df) 
    }
  df[is.na(df[[colname]]) | df[[colname]] == "." | suppressWarnings(as.numeric(df[[colname]])) <= threshold, ]
}
filter_AF_more <- function(df, colname, threshold = 0.05) {
  if (!colname %in% names(df)) {
    warning(paste("Column", colname, "not found in the data frame. Function skipped."))
    return(df) 
  }
  df[is.na(df[[colname]]) | df[[colname]] == "." | suppressWarnings(as.numeric(df[[colname]])) >= threshold, ]
}

# Stepwise Filtering
Variants1 <- subset(Variants, !grepl("benign", CLNSIG, ignore.case = TRUE))
Variants1 <- filter_AF_less(Variants1, "X1000g2015aug_all")
Variants1 <- filter_AF_less(Variants1, "gnomad40_exome_AF")
Variants1 <- filter_AF_less(Variants1, "gnomad40_genome_AF")
Variants1 <- filter_AF_less(Variants1, "Kaviar_AF")
Variants1 <- filter_AF_less(Variants1, "esp6500siv2_all")
Variants1 <- filter_AF_less(Variants1, "GME_AF", threshold = 0.1)
Variants1 <- filter_AF_less(Variants1, "HRC_AF")
Variants1 <- filter_AF_less(Variants1, "Iranome_AF", threshold = 0.1)

# Functional Annotation Filters
Variants1 <- subset(Variants1,
  grepl("exonic|splicing", Func_refGeneWithVer, ignore.case = TRUE) |
  grepl("exonic|splicing", Func_ensGene, ignore.case = TRUE) |
  grepl("exonic|splicing", Func_knownGene, ignore.case = TRUE)
)

Variants1 <- subset(Variants1,
  ExonicFunc_refGeneWithVer != "synonymous SNV" &
  ExonicFunc_ensGene != "synonymous SNV" &
  ExonicFunc_knownGene != "synonymous SNV"
)

# Interpretation & Prediction Score Filters
Variants1 <- subset(Variants1, !grepl("benign", InterVar_automated, ignore.case = TRUE))

Variants1 <- filter_AF_more(Variants1, "DANN_score", threshold = 0.75)
Variants1 <- filter_AF_more(Variants1, "REVEL", threshold = 0.25)

# Generate output file name
input_base <- tools::file_path_sans_ext(basename(input_file))
output_dir <- dirname(input_file)
output_file <- file.path(output_dir, paste0(input_base, "_filtered.csv"))
# Write Output
write.csv(Variants1, file = output_file, quote = FALSE, row.names = FALSE)

cat("✅ Filtering done!\n📄 Output file:", output_file, "\n")
