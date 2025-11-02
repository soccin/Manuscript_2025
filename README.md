# Single-Cell Whole-Genome Sequencing Analysis Pipeline

Analysis pipeline for detecting and characterizing discordant read pairs from single-cell DNA sequencing data to identify structural variants and chromosomal rearrangements.

## Overview

This pipeline processes aligned single-cell whole-genome sequencing data (BAM files) to extract properly paired reads and map them to fixed-width genomic bins, facilitating genome-wide analysis of read pair connectivity patterns and structural variant detection.

## Pipeline Steps

### Step 1: Paired-End Read Processing and Genomic Bin Assignment

**Script:** `bin/bam2bedpe.sh`

Extracts high-quality paired-end reads from BAM files and maps them to 20 kb genomic bins across the hg38 reference genome.

**Usage:**
```bash
./bin/bam2bedpe.sh <input.bam>
```

**Output:**
- `out/20k/<dir>/<sample>.pe.bed.gz` - Filtered paired-end reads in BEDPE format
- `out/20k/<dir>/<sample>.pe.map.gz` - Read pairs mapped to genomic bins

**What it does:**
1. Filters reads for quality (removes unmapped, secondary, duplicate, and QC-failed reads)
2. Converts to BEDPE format (both read ends on one line)
3. Maps each read end to a 20 kb genomic bin
4. Retains only read pairs where both ends map to defined bins

## Requirements

### Software Dependencies

- **SAMtools** (v1.x+) - BAM file processing
- **BEDTools** (v2.x+) - BEDPE format conversion
- **R** (v4.x+) - Statistical analysis
  - tidyverse
  - tidygenomics

### Data Requirements

- Aligned BAM files from single-cell whole-genome sequencing
- Genomic bins file: `data/raw/bins/hg38_20k_gz_enc_bins.bed`

## Directory Structure

```
.
├── bin/
│   ├── bam2bedpe.sh          # Main pipeline script
│   └── src/
│       └── map2bins.R        # R script for genomic bin assignment
├── data/
│   └── raw/
│       └── bins/             # Genomic bin definitions
├── out/                      # Output directory (auto-created)
│   └── 20k/                  # 20kb bin analysis results
├── tests/                    # Unit tests
│   └── UNITTEST02.sh         # Test for bam2bedpe pipeline
└── manu/
    └── methods/              # Methods documentation for manuscript
```

## Testing

Run unit tests to verify pipeline functionality:

```bash
./tests/UNITTEST02.sh         # Quick test
FULL=Yes ./tests/UNITTEST02.sh  # Full test suite
```

## Documentation

Detailed methods documentation suitable for publication is available in:
- `manu/methods/STEP_01_v1.md` - Complete technical methods
- `manu/methods/STEP_01_SUMMARY_v1.md` - Condensed narrative version

## Example

```bash
# Process a single-cell BAM file
./bin/bam2bedpe.sh data/raw/bams/sample001.bam

# Output will be created in:
# out/20k/data/raw/bams/sample001.pe.bed.gz
# out/20k/data/raw/bams/sample001.pe.map.gz
```

## Citation

If you use this pipeline, please cite:

[Citation information to be added upon publication]

## Contact

[Contact information]

## License

[License information]
