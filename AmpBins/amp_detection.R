#' Annotate Amplification Block Cluster
#'
#' Combines genomic coordinates, cytoband, genes, and sample metadata for a
#' cluster of overlapping amplification blocks.
#'
#' @param block_cluster Tibble containing amplification blocks from one cluster
#'
#' @return Tibble with cluster summary including genomic location, cytobands,
#'   genes, and sample information
annotate_block <- function(block_cluster) {
  start <- min(block_cluster$start)
  end <- max(block_cluster$end)

  # Extract genomic coordinates and cytoband information
  ginfo <- lpsA_chromInfo |> slice(start:end)
  cyto <- ginfo |>
    summarize(across(c(arm, band), ~paste(unique(.x), collapse = ";")))

  gStartMb <- round(min(ginfo$bin.start) / 1e6, 2)
  gEndMb <- round(max(ginfo$bin.end) / 1e6, 2)

  # Extract sample metadata
  sample_metadata <- block_cluster |>
    left_join(select(lpsA_Y, cellID, subsample, bioID, subtype, clone.type),by=join_by(cellID)) |>
    select(subtype, clone.type, bioID, subsample) |>
    summarize(across(everything(), ~paste0(unique(.x), collapse = ";")))

  # Extract genes in this region
  genes <- lpsA_gene.index |>
    select(bin.id, Gene = hgnc.symbol) |>
    slice(start:end) |>
    summarize(Genes = paste(unique(Gene), collapse = ";"))

  # Combine all annotations
  tibble(
    cluster = block_cluster$cluster[1],
    NCells = block_cluster$n_cluster[1],
    NUBlocks = block_cluster$NUniq[1],
    start = start,
    end = end
  ) |>
    bind_cols(cyto) |>
    mutate(gStart = gStartMb, gEnd = gEndMb) |>
    bind_cols(sample_metadata) |>
    bind_cols(genes)
}

#' Detect Amplifications in Copy Number Data
#'
#' Identifies contiguous regions where copy number values exceed a threshold.
#'
#' @param x Numeric vector of copy numbers to search
#' @param cutoff Copy number threshold value
#' @param min_length Minimum consecutive bins needed (default 3)
#'
#' @return List containing:
#'   - islands: Data frame with start and end indices of each amplification
#'   - runs: RLE object for downstream analysis
#'
#' @examples
#' cn <- c(1, 5, 6, 7, 2, 8, 9, 10, 11, 3, 12, 13, 14)
#' detect_amplifications(cn, cutoff = 5)
#' # $islands
#' #   start end
#' # 1     2   4
#' # 2     6   9
#' # 3    11  13
detect_amplifications <- function(x, cutoff, min_length) {
  above_cutoff <- x > cutoff

  # Use run-length encoding to find consecutive TRUE values
  runs <- rle(above_cutoff)

  # Identify runs that exceed both cutoff and minimum length
  valid_runs <- runs$lengths >= min_length & runs$values

  if (!any(valid_runs)) {
    return(list(
      islands = data.frame(start = integer(0), end = integer(0)),
      runs = runs
    ))
  }

  # Calculate start and end positions for each run
  end_positions <- cumsum(runs$lengths)
  start_positions <- c(1, end_positions[-length(end_positions)] + 1)

  amplifications <- data.frame(
    start = start_positions[valid_runs],
    end = end_positions[valid_runs]
  )

  list(islands = amplifications, runs = runs)
}

#' Summarize Amplifications with Statistics
#'
#' Detects amplified regions and calculates summary statistics for each.
#'
#' @param Xi Single-column tibble containing copy number values
#' @param cutoff Copy number threshold for identifying amplifications
#' @param min_length Minimum consecutive bins needed (default 3)
#'
#' @return Tibble with amplification boundaries and statistics (start, end, N,
#'   mean, sd, min, max, island number, cellID)
summarize_amplifications <- function(Xi, cutoff, min_length = 3) {
  cellID <- colnames(Xi)[1]
  copy_numbers <- Xi[[1]]

  # Detect amplified regions
  amp_result <- detect_amplifications(copy_numbers, cutoff, min_length)

  # Calculate statistics for each rle block
  runs <- amp_result$runs
  block_id <- rep(seq_along(runs$lengths), runs$lengths)

  amp_stats <- tibble(tcn = copy_numbers, block = block_id) |>
    group_by(block) |>
    summarize(
      N = n(),
      mean = mean(tcn),
      sd = sd(tcn),
      min = min(tcn),
      max = max(tcn)
    ) |>
    mutate(is_amplified = runs$values) |>
    filter(is_amplified) |>
    mutate(island = row_number())

  # Combine with amplification boundary information
  bind_cols(amp_result$islands, amp_stats) |>
    select(-block, -is_amplified) |>
    mutate(cellID = cellID)
}
