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
Variants <- read.delim(file = input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("✅ Reading input file finished!\n")
# A function to read allele frequencies safely
parse_af <- function(x) {
  x_clean <- trimws(x)
  x_clean[x_clean %in% c("", ".")] <- NA              # explicit missings
  v <- strsplit(x_clean, "[,;|]")                     # multi-allelic or sub-pop
  sapply(v, function(z) {
    z_num <- suppressWarnings(as.numeric(z))
    if (all(is.na(z_num))) NA else max(z_num, na.rm = TRUE)
  })
}

# Helper to safely convert to numeric and filter based on threshold
filter_AF_less <- function(df, colname, threshold = 0.05) {
  if (!colname %in% names(df)) {
    warning(sprintf("Column %s not found – skipped.", colname)); return(df)
  }
  af <- parse_af(df[[colname]])
  df[is.na(af) | af <= threshold, , drop = FALSE]
}

filter_AF_more <- function(df, colname, threshold = 0.05) {
  if (!colname %in% names(df)) {
    warning(sprintf("Column %s not found – skipped.", colname)); return(df)
  }
  af <- parse_af(df[[colname]])
  df[!is.na(af) & af >= threshold, , drop = FALSE]
}

# Stepwise Filtering
Variants1 <- subset(Variants, !grepl("benign", CLNSIG, ignore.case = TRUE))
Variants1 <- filter_AF_less(Variants1, "gvs_all_af")
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

Variants1 <- filter_AF_more(Variants1, "DANN_score", threshold = 0.8)
Variants1 <- filter_AF_more(Variants1, "REVEL", threshold = 0.3)

# Generate output file name
input_base <- tools::file_path_sans_ext(basename(input_file))
output_dir <- dirname(input_file)
output_file <- file.path(output_dir, paste0(input_base, "_filtered.tsv"))
# Write Output
write.table(Variants1, file = output_file, quote = FALSE, row.names = FALSE, sep = "\t")

cat("✅ Filtering done!\n📄 Output file:", output_file, "\n")
