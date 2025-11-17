# ==============================================================================
# Amplification Block Analysis
# ==============================================================================
# Detect and analyze focal amplifications in single-cell copy number data,
# excluding chromosome 12q amplifications

source("amp_detection.R")
require(tidyverse)
require(furrr)

# Load reference data
data(lpsA_X)
data(lpsA_Y)
data(lpsA_chromInfo)
data(lpsA_gene.index)

# Clone type counts for reference:
# EDC: 224, Adv.Tumor: 5188, Adv.Lipoma: 436
# HC: 1859, Stroma: 7308, Immune: 2158

# ==============================================================================
# Select Cell Subset
# ==============================================================================

ALL = FALSE
if (ALL) {
  cells <- lpsA_Y$cellID
} else {
  cells <- lpsA_Y |>
    filter(clone.type %in% c("EDC", "Adv.Tumor"),
           subtype != "Lipoma") |>
    pull(cellID)
}

X <- lpsA_X |>
  select(bin.id, all_of(cells))

# ==============================================================================
# Detect Amplifications
# ==============================================================================

# Use parallel processing for speed
plan(multisession, workers = 16)

# Detect amplifications across all cells
amp_blocks <- X |>
  select(-bin.id) |>
  future_imap_dfr(
    ~summarize_amplifications(tibble(!!.y := .x), cutoff = 10, min_length = 3),
    .progress = TRUE
  )

# Filter out amplifications on chromosome arm 12q
bins_12q <- lpsA_chromInfo |>
  filter(arm == "12q") |>
  pull(bin.id)

amp_blocks_not_12q <- amp_blocks |>
  filter(!(start %in% bins_12q & end %in% bins_12q))

# ==============================================================================
# Summarize Per-Cell Statistics
# ==============================================================================

cell_stats <- amp_blocks_not_12q |>
  group_by(cellID) |>
  summarize(n_blocks = n()) |>
  right_join(tibble(cellID = cells)) |>
  mutate(n_blocks = replace_na(n_blocks, 0)) |>
  arrange(n_blocks) |>
  left_join(select(lpsA_Y, cellID, sample.type, clone.type, subtype))

# ==============================================================================
# Plot Distribution
# ==============================================================================

plot_base <- cell_stats |>
  ggplot(aes(n_blocks, fill = paste0(subtype, "|", clone.type))) +
  geom_histogram(binwidth = 1, alpha = 0.5, color = "grey35")

plot_faceted <- plot_base +
  facet_grid(clone.type ~ subtype, scale = "free_y") +
  theme_light(16) +
  scale_fill_brewer(palette = "Paired", name = "SubType|Clone") +
  guides(fill = guide_legend(reverse = TRUE))

pdf(file = "amplificationBlocks_TCN_11_Length_3_v2.pdf", width = 11, height = 8.5)
print(plot_faceted)
dev.off()

# ==============================================================================
# Cluster Overlapping Amplifications
# ==============================================================================

# Identify overlapping amplification blocks across cells
# Algorithm: increment cluster ID when a block starts after the max end of previous blocks
amp_cluster <- amp_blocks_not_12q |>
  arrange(start, end) |>
  mutate(cluster = cumsum(start > lag(cummax(end), default = -1L))) |>
  group_by(cluster) |>
  mutate(n_cluster = n()) |>
  ungroup() |>
  group_by(start, end) |>
  mutate(NUniq = n()) |>
  ungroup()

# ==============================================================================
# Identify Sporadic Amplifications
# ==============================================================================
# Extract and annotate low-frequency amplifications (≤5 cells) to identify
# potential private or early clonal events distinct from widespread alterations

sporadic_blocks <- amp_cluster %>% filter(n_cluster <= 5) %>% split(.$cluster)
sporadic_clusters <- map_dfr(sporadic_blocks, annotate_block)
# Create a list of data frames for multiple sheets
output_data <- list(
  SporadicClusters = sporadic_clusters,
  ColumnDescriptions = read_csv("sporadic_clusters_column_descriptions.csv")
)

write_xlsx(output_data, "sporadicClusters_TCN_11_Len_3_Size_5_v1.xlsx")

# Sporadic Clusters Output Columns:
# ----------------------------------
# cluster     - Unique identifier for each genomic region with overlapping amplifications
# NCells      - Number of cells carrying an amplification in this region
# NUBlocks    - Number of distinct amplification boundaries within this cluster
#               (same region may be amplified with slightly different start/end points)
# start       - Genomic bin index marking the earliest amplification start
# end         - Genomic bin index marking the latest amplification end
# arm         - Chromosome arm location (1p, 1q, etc.); semicolon-separated if spanning multiple
# band        - Cytoband annotation (e.g., 1p36.21); semicolon-separated for multi-band regions
# gStart      - Genomic start coordinate in megabases (Mb)
# gEnd        - Genomic end coordinate in megabases (Mb)
# subtype     - Tumor histological subtype(s) carrying this amplification
# clone.type  - Clonal classification(s) of cells with this amplification
# bioID       - Biological sample identifier(s)
# subsample   - Specific subsample identifier(s)
# Genes       - Gene symbols within the amplified region; semicolon-separated list
#
# Note: Metadata columns (subtype, clone.type, bioID, subsample) are semicolon-
#       separated when multiple cells from different samples share the amplification
