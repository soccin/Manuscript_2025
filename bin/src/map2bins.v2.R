usage = function() {
  cat("
Usage: Rscript map2bins.R <bins.bed> <mapped_reads.bed.gz>

Description:
  Map paired-end reads to genomic bins and filter for valid pairs.

Arguments:
  bins.bed           - Tab-delimited file with genomic bins (chr, start, end)
  mapped_reads.bed.gz - Tab-delimited BED file with paired-end mappings
                        (chr1, start1, end1, chr2, start2, end2)

Output:
  Creates a .map.gz file containing mappings with BIN1 and BIN2 columns.

")
  quit(status = 1)
}

argv = commandArgs(trailing = TRUE)

if (length(argv) != 2) {
  usage()
}

BINFILE = argv[1]
MAPFILE = argv[2]

suppressPackageStartupMessages({
  require(tidygenomics)
  require(tidyverse)
})

# Suppress lifecycle warnings from tidygenomics
options(lifecycle_verbosity = "quiet")

# Read genomic bins and assign bin IDs
bins = read_tsv(BINFILE, col_names = FALSE, show_col_types=F) |>
  mutate(X2 = X2 + 1) |>
  mutate(BIN = sprintf("B%05d", row_number()))

# Read paired-end mapping data
map = read_tsv(MAPFILE, col_names = FALSE, show_col_types=F)

# Check for empty input
if (nrow(map) == 0) {
  cat("\nERROR: No PE READS passed for", MAPFILE, "\n\n")
  quit(status = 1)
}

# Create end coordinates and assign mapping IDs
map = map |>
  mutate(X2=X2+1) |>
  mutate(X3 = X2, X5 = X6) |>
  mutate(MID = sprintf("M%05d", row_number()))

# Find bin overlaps for both ends of paired reads
bin1 = genome_intersect(map, bins, by = c(X1 = "X1", X2 = "X2", X3 = "X3")) |>
  select(MID, BIN1 = BIN)

bin2 = genome_intersect(map, bins, by = c(X4 = "X1", X5 = "X2", X6 = "X3")) |>
  select(MID, BIN2 = BIN)

# Join bin assignments and filter for valid pairs
map = map |>
  left_join(bin1, by = "MID") |>
  left_join(bin2, by = "MID") |>
  filter(!is.na(BIN1) & !is.na(BIN2))

# Write output
output_file = gsub(".bed.gz", ".map.gz", MAPFILE)
write_tsv(map, output_file)

