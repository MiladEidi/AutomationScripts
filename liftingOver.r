library(rtracklayer)
library(GenomicRanges)

# File paths
input_file <- "/lustre06/project/6003434/miladad/hg38_gnomad41_genome.txt.gz"
output_mapped <- "/lustre06/project/6003434/miladad/hg19_gnomad41_genome.txt"
output_unmapped <- "/lustre06/project/6003434/miladad/hg19_gnomad41_genome_Unmapped.txt"

# Open the gzipped file and read the first line (header)
gz_con <- gzfile(input_file, "rt")  # Open the file in text mode
header <- readLines(gz_con, n = 1) # Read the first line

header_cols <- unlist(strsplit(header, "\t"))
updated_header <- paste(header_cols, collapse = "\t")

# Write the updated header to the output files
writeLines(updated_header, output_mapped)
writeLines(updated_header, output_unmapped)

# Chain file
chain <- import.chain("/lustre06/project/6003434/miladad/hg38ToHg19.over.chain")

# Process the input in chunks
chunk_size <- 100000
readLines(gz_con, n = 1) # Skip the header line in subsequent reads

while (TRUE) {
  # Read a chunk of data
  table_chunk <- read.table(gz_con, nrows = chunk_size, header = FALSE, sep = "\t")
  
  if (nrow(table_chunk) == 0) break # Stop if no rows are read
  
  # Create GRanges object
  coords <- GRanges(seqnames = table_chunk[, 1], ranges = IRanges(start = table_chunk[, 2], end = table_chunk[, 3]))
  mcols(coords) <- table_chunk[, -c(1:3)]
  
  # Perform liftover
  lifted <- liftOver(coords, chain)
  mapped <- lengths(lifted) > 0
  
  # Write mapped positions (excluding third and fourth columns)
  lifted_flat <- as.data.frame(unlist(lifted[mapped]))[, -c(4, 5)]
  write.table(lifted_flat, file = output_mapped, sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Write unmapped positions (excluding third and fourth columns)
  unmapped_coords <- as.data.frame(coords[!mapped])[, -c(4, 5)]
  write.table(unmapped_coords, file = output_unmapped, sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Close the file connection
close(gz_con)

