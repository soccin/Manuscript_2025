## Detection of Discordant Read Pairs from Single-Cell Whole-Genome Sequencing

To identify structural variants and chromosomal rearrangements in single-cell DNA sequencing data, we developed a pipeline to detect and characterize discordant read pairs—defined as paired-end reads mapping either to different chromosomes or separated by distances substantially exceeding the expected insert size (~20 kb).

**Read Pair Extraction and Quality Filtering**

For each single-cell BAM file containing aligned whole-genome sequencing reads, we first extracted high-quality properly paired reads using SAMtools (v1.x) and BEDtools (v2.x). BAM files were sorted by read name to ensure mate pairs were adjacent, then filtered to retain only properly paired reads while excluding unmapped reads, reads with unmapped mates, secondary alignments, quality control failures, and PCR duplicates (SAMtools flags: -f 1 -F 3084). Filtered read pairs were converted to BEDPE format using BEDtools bamtobed, producing six-column records containing the genomic coordinates of both read ends (chr1, start1, end1, chr2, start2, end2).

**Genomic Binning and Read Pair Assignment**

To facilitate genome-wide analysis of read pair connectivity patterns, we mapped paired-end reads to fixed-width genomic bins. The hg38 reference genome was divided into non-overlapping 20 kb bins spanning all autosomes and sex chromosomes. Each read pair was assigned to two bins (BIN1 and BIN2) based on the genomic overlap of each read end with the bin coordinates. Bin assignment was performed using the tidygenomics R package, which implements efficient genomic interval intersection operations. Read pairs were retained for downstream analysis only if both read ends mapped unambiguously to defined genomic bins. Each bin was assigned a unique identifier (B00001, B00002, etc.) for subsequent connectivity analysis.

This binning strategy enabled systematic quantification of read pair connectivity between genomic regions, facilitating the identification of recurrent structural variants and the reconstruction of chromosomal rearrangement patterns across the single-cell population.
