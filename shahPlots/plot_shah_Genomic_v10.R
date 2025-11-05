# Structural Variant Link Analysis
# Extracted from interactive session: work01.R

# Load required libraries and functions
source("load_shah.R")
source("tools.R")
source("plotTools.R")
require(tidyverse)

# Configuration
SAMPLE_ID <- "WD0539P"
TCN_FILTER <- 10  # Total copy number filter threshold
PE_FILTER <- 3    # Paired-end reads filter threshold
NUM_CELLS <- 10

#
data(lpsA_chromInfo)

# Load sample data (all clones)
sample_data <- load_shah_allClones(SAMPLE_ID)

# Explore data structure
cat("\nSample data structure:\n")
cat("  Y (metadata):   ", nrow(sample_data$Y), "cells\n")
cat("  maps (SVs):     ", nrow(sample_data$maps), "links\n\n")

# Show clone distribution
cat("Clone distribution:\n")
cloneCounts=sample_data$Y |>
  count(final_cluster) |>
  arrange(desc(n))

clones=cloneCounts %>% filter(n>=NUM_CELLS) %>% pull(final_cluster)

clone=clones[1]

pc=list()
linkStats=list()

for(clone in clones) {
  cat("Clone =",clone,"...")
  clone_subset=subset_cells(sample_data,final_cluster==clone)
  clone_links=get_links(clone_subset,TCN_FILTER,PE_FILTER)
  link_summary=summarize_links(clone_links,SAMPLE_ID,TCN_FILTER,PE_FILTER,clone)
  p1=plot_genomic_links(clone_links,SAMPLE_ID,clone)+geom_hline(yintercept=c(TCN_FILTER),color="darkcyan",alpha=.5)
  linkStats[[clone]]=link_summary
  pc[[clone]]=p1
  cat("\n")
}

pfile=cc("shahGenomic",SAMPLE_ID,"10CellClones","TCNFilt",TCN_FILTER,"PEFilter",PE_FILTER,"v10.pdf")
pdf(file=pfile,width=14,height=8.5)
print(pc)
dev.off()