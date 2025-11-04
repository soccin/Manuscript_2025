# Cache database for memoised functions
cache_db <- cachem::cache_disk("./cache")

# Configuration
MAPDIR <- "data/raw/maps/v2023"

#' Select standard Y metadata columns
#'
#' @param data Data frame to select columns from
#' @return Data frame with selected columns
select_y_columns <- function(data) {
  data |>
    select(
      cellID:subtype,
      barcode,
      ConsensusC,
      final_cluster,
      spatial.extent,
      clone.type,
      Branch,
      singerNo
    )
}


#' Validate sample exists and show available samples if not
#'
#' @param Y Filtered metadata data frame
#' @param SampleID Sample identifier that was searched for
validate_sample <- function(Y, SampleID) {
  if (nrow(Y) == 0) {
    data(lpsA_Y)
    valid_samples <- lpsA_Y |> distinct(sampleID) |> pull()
    cat("\n\tNo Samples matching", SampleID, "\n")
    cat("\tValid samples:", valid_samples, "\n\n")
    quit()
  }
}


#' Load copy number data for specific cells
#'
#' @param cell_ids Vector of cell IDs to load
#' @return Data frame with bin.id, cellID, and tcn columns in long format
load_copy_number <- function(cell_ids) {
  data(lpsA_X)
  lpsA_X |>
    select(bin.id, all_of(cell_ids)) |>
    gather(cellID, tcn, -bin.id)
}


#' Load and process structural variant maps
#'
#' @param bioID Biological sample ID for file lookup
#' @param cell_ids Vector of cell IDs to filter
#' @return Processed maps data frame
load_sv_maps <- function(bioID, cell_ids) {
  MAPFILE <- file.path(MAPDIR, cc("maps", bioID, ".rds"))

  readRDS(MAPFILE) |>
    filter(cellID %in% cell_ids) |>
    filter(BIN1 != BIN2) |>
    mutate(
      BIN1 = gsub("B", "", BIN1) |> as.numeric(),
      BIN2 = gsub("B", "", BIN2) |> as.numeric()
    ) |>
    filter(abs(BIN1 - BIN2) > 1) |>
    group_by(CHR1, BIN1, CHR2, BIN2, N, Nx, cellID) |>
    summarize(n = sum(n), .groups = "drop")
}


#' Load Shah data for a specific sample and clone type
#'
#' @param SampleID Sample identifier to load
#' @param CloneType Clone type to filter (e.g., "Adv", "NonAdv")
#' @return List with three elements: X (copy number data), Y (cell metadata),
#'   maps (structural variant maps)
.load_shah <- function(SampleID, CloneType) {

  # Load and filter cell metadata
  data(lpsA_Y)
  Y <- lpsA_Y |>
    filter(sampleID == SampleID & clone.type == CloneType) |>
    select_y_columns()

  validate_sample(Y, SampleID)

  # Load copy number and structural variant data
  X <- load_copy_number(Y$cellID)
  maps <- load_sv_maps(Y$bioID[1], Y$cellID)

  list(X = X, Y = Y, maps = maps)
}

# Memoised version for performance
load_shah <- memoise::memoise(.load_shah, cache = cache_db)


#' Load Shah data for all clones in a sample
#'
#' @param SampleID Sample identifier to load
#' @return List with three elements: X (copy number data), Y (cell metadata),
#'   maps (structural variant maps with clone annotations)
.load_shah_allClones <- function(SampleID) {

  # Load and filter cell metadata (all clones with consensus clustering)
  data(lpsA_Y)
  Y <- lpsA_Y |>
    filter(sampleID == SampleID & !is.na(ConsensusC)) |>
    select_y_columns()

  validate_sample(Y, SampleID)

  # Load copy number and structural variant data
  X <- load_copy_number(Y$cellID)
  maps <- load_sv_maps(Y$bioID[1], Y$cellID) |>
    left_join(lpsA_Y |> select(cellID, ConsensusC, clone.type), by = "cellID")

  list(X = X, Y = Y, maps = maps)
}

# Memoised version for performance
load_shah_allClones <- memoise::memoise(.load_shah_allClones, cache = cache_db)
