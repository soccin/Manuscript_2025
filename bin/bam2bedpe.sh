#!/bin/bash
#
# bam2bedpe.sh - Convert BAM file to BEDPE format and map to genomic bins
#
# Usage: bam2bedpe.sh <input.bam>
#
# Description:
#   Processes a BAM alignment file to extract properly paired reads,
#   converts them to BEDPE format (both ends of paired reads), and
#   maps the read pairs to predefined genomic bins.
#
# Output:
#   - out/<dir>/<basename>.pe.bed.gz - Paired-end reads in BEDPE format
#   - out/<dir>/<basename>.pe.map.gz - Reads mapped to bins (from map2bins.R)
#

usage() {
  cat << EOF

Usage: bam2bedpe.sh <input.bam>

Description:
  Extract paired-end reads from BAM file and map to genomic bins.

Arguments:
  input.bam - Aligned BAM file to process

Output:
  Creates BEDPE file and bin-mapped file in out/<dir>/ directory.

EOF
  exit 1
}

# Check for required argument
if [ $# -ne 1 ]; then
  usage
fi

# Get script directory for finding helper scripts
SDIR=$(dirname "$(readlink -f "$0")")

# Load required modules
module load samtools
module load bedtools

# Define genomic bins file (20kb bins for hg38)
BINS=data/raw/bins/hg38_20k_gz_enc_bins.bed
BTAG=$(basename ${BINS/_gz_*/})

# Parse input BAM file path
BAM=$1
BASE=$(basename ${BAM/.bam/})
DIR=$(dirname $BAM)

# Create output directory
ODIR=out/20k/$DIR
mkdir -p $ODIR

# Process BAM file to extract paired-end reads
# Pipeline steps:
#   1. Sort by read name (-n) to pair mates together
#   2. Filter reads:
#      -F 3084: Exclude unmapped, secondary, QC-fail, and duplicate reads
#               (3084 = 0x4 + 0x8 + 0x100 + 0x400 + 0x800)
#      -f 1:    Keep only paired reads
#   3. Convert to BEDPE format (both ends of pair on one line)
#   4. Replace spaces with tabs for clean TSV format
#   5. Compress output
samtools sort -n $BAM -@ 16 -m 1g \
    | samtools view -F 3084 -f 1 -b \
    | bedtools bamtobed -bedpe -i - \
    | tr ' ' '\t' \
    | gzip -9 - \
    > $ODIR/${BASE}.pe.bed.gz

# Map paired-end reads to genomic bins
Rscript $SDIR/src/map2bins.R $BINS $ODIR/${BASE}.pe.bed.gz

