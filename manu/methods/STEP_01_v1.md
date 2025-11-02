# Methods Draft: Paired-End Read Processing and Genomic Bin Assignment

## Processing of Paired-End Sequencing Data

### Genome Binning Strategy

The human reference genome (hg38) was partitioned into non-overlapping 20 kb bins, generating a total of [N] genomic bins spanning all autosomal and sex chromosomes. This bin size was selected to balance spatial resolution with sufficient coverage per bin for robust statistical analysis. Each bin was assigned a unique identifier for downstream tracking of chromatin interactions.

### Read Quality Filtering and Processing

Aligned paired-end sequencing reads (BAM format) were processed through a multi-step quality control and formatting pipeline. First, reads were sorted by read name using SAMtools (v[VERSION]) to ensure proper pairing of read mates. Stringent quality filters were then applied to retain only high-confidence alignments: reads were excluded if they were unmapped, had an unmapped mate, represented secondary alignments, failed quality control checks, or were marked as PCR duplicates (SAMtools flags: -F 3084 -f 1). This filtering strategy ensured that only properly paired, primary alignments from unique molecules were retained for analysis.

The filtered paired-end reads were converted to BEDPE (Browser Extensible Data Paired-End) format using BEDTools (v[VERSION]), which represents both ends of each read pair on a single line with their respective genomic coordinates (chromosome, start, end positions for both mates). This format facilitates efficient downstream processing of read-pair relationships.

### Assignment of Read Pairs to Genomic Bins

Each paired-end read was assigned to genomic bins based on the mapping positions of both read ends. For each read pair, we independently determined which 20 kb bin overlapped with the first read end (Read 1) and which bin overlapped with the second read end (Read 2), using genome interval intersection operations implemented in the tidygenomics R package (v[VERSION]).

Read pairs were retained for downstream analysis only if both read ends successfully mapped to defined genomic bins. Read pairs where one or both ends fell outside bin boundaries (e.g., at chromosome ends or in unmappable regions) were excluded. Each valid read pair was annotated with two bin identifiers (BIN1 and BIN2), representing the genomic bins contacted by the two ends of the sequenced DNA fragment, thereby creating a discrete representation of potential chromatin interactions.

### Quality Control

Samples that yielded no valid paired-end reads passing all filtering criteria were flagged and excluded from downstream analysis. For each sample, we recorded the total number of input read pairs, the number passing quality filters, and the number successfully assigned to genomic bins to assess data quality and processing efficiency.

### Software and Computational Implementation

Read processing was performed using SAMtools v[X.X.X] with 16 parallel threads and 1 GB memory per thread. Genomic interval operations used BEDTools v[X.X.X]. Bin assignment and statistical analysis were conducted in R v[X.X.X] using the tidyverse (v[X.X.X]) and tidygenomics (v[X.X.X]) packages. All analysis scripts are available at [REPOSITORY URL].

---

## Notes for Completion:

1. **Fill in software versions**: Check actual versions used in your environment
   - SAMtools version: `samtools --version`
   - BEDTools version: `bedtools --version`
   - R version: `R --version`
   - R package versions: `packageVersion("tidyverse")`, etc.

2. **Fill in [N] bins**: Calculate total number of 20kb bins in hg38
   - Approximately ~154,000 bins for hg38 (3.1 Gb / 20kb)

3. **Add context**: This section assumes this is the first step in a larger analysis pipeline. You may need to add:
   - What type of data this is (Hi-C, ChIA-PET, proximity ligation, etc.)
   - How the original BAM files were generated (alignment parameters, reference)
   - What happens to this data in subsequent analysis steps

4. **Repository**: Add your GitHub/GitLab repository URL when publishing

5. **Consider adding**: Depending on journal and reviewers:
   - Specific ENCODE or reference for bin definitions if using standard bins
   - Rationale for 20kb bin size (cite similar studies)
   - Statistics on typical filtering rates in your data
