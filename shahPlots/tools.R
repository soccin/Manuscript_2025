#' Subset cells from Shah data object
#'
#' Filters all components (X, Y, maps) of a Shah data object based on
#' a predicate expression applied to the Y metadata table.
#'
#' @param oo Shah data object (list with X, Y, maps elements)
#' @param predicate Filtering expression (unquoted, applied to Y table)
#' @return Modified Shah data object with filtered cells, includes
#'   `subset_filter` element containing the filter expression as text
#'
#' @examples
#' subset_cells(oo, clone.type == "EDC")
#' subset_cells(oo, final_cluster == "WD0539.X5")
#' subset_cells(oo, clone.type == "EDC" & spatial.extent > 0.5)
subset_cells <- function(oo, predicate) {

  # Get cells matching the predicate
  cells <- oo$Y |>
    filter({{ predicate }}) |>
    pull(cellID)

  # Filter all data components to matching cells
  oo$X <- oo$X |> filter(cellID %in% cells)
  oo$Y <- oo$Y |> filter(cellID %in% cells)
  oo$maps <- oo$maps |> filter(cellID %in% cells)

  # Store filter expression for reference
  oo$subset_filter <- rlang::quo_text(rlang::enquo(predicate))

  oo
}

#' Extract and filter structural variant links
#'
#' Processes structural variant maps to identify significant links
#' based on copy number and paired-end read thresholds.
#'
#' @param os Shah data object (output from load_shah or subset_cells)
#' @param TCN.FILT Total copy number filter threshold
#' @param PE.FILT Paired-end read count filter threshold
#' @return List with elements:
#'   - links: Filtered structural variant links with type annotations
#'   - cnv: Median copy number by bin
#'   - cellNumScaling: Normalization factor (scales to 100 cells)
#'   - tcnFilterBins: Bins passing copy number threshold
get_links <- function(os, TCN.FILT, PE.FILT) {

  # Normalize for cell count by scaling to 100 fixed cells
  cellNumScaling <- 100 / nrow(os$Y)

  # Aggregate and normalize structural variant maps
  maps <- os$maps |>
    mutate(NC = cellNumScaling * n / (N / 1e6)) |>
    group_by(CHR1, BIN1, CHR2, BIN2) |>
    summarize(n = sum(NC), .groups = "drop") |>
    mutate(Type = case_when(
      CHR1 == "chr12" & CHR2 == "chr12" ~ "12->12",
      CHR1 == "chr12" | CHR2 == "chr12" ~ "12->A",
      CHR1 != CHR2 ~ "A->B",
      TRUE ~ "A->A[!12]"
    ))

  # Calculate median copy number per bin
  cnv <- os$X |>
    mutate(bin.id = as.numeric(bin.id)) |>
    group_by(bin.id) |>
    summarize(med.tcn = median(tcn), .groups = "drop")

  # Identify bins with copy number above threshold
  tcnFilterBins <- cnv |>
    filter(med.tcn >= TCN.FILT) |>
    pull(bin.id)

  # Filter links by read count and copy number thresholds
  links <- maps |>
    filter(n >= PE.FILT) |>
    filter(BIN1 %in% tcnFilterBins & BIN2 %in% tcnFilterBins)

  list(
    links = links,
    cnv = cnv,
    nCells = nrow(os$Y),
    cellNumScaling = cellNumScaling,
    tcnFilterBins = tcnFilterBins
  )
}

#' Get bins of interest for genomic plotting
#'
#' Identifies genomic bins for visualization based on structural variant
#' links and copy number alterations.
#'
#' @param sid Sample ID
#' @param size Maximum cumulative width in bases (default: 100Mb)
#' @param ct Clone type
#' @return Vector of bin IDs representing regions of interest
get_bins_of_interest <- function(sid, size = 100e6, ct) {

  ob <- get_links(sid, ct)

  # Extract bins involved in inter-chromosomal links
  link_bins <- ob$links |>
    filter(CHR1 != CHR2) |>
    select(BIN1, BIN2) |>
    unlist() |>
    unname()

  # Find high copy number bands outside chr12q that are linked
  band_rois <- ob$cnv |>
    left_join(lpsA_chromInfo, by = "bin.id") |>
    arrange(desc(med.tcn)) |>
    filter(arm != "12q") |>
    distinct(band, .keep_all = TRUE) |>
    mutate(linked = bin.id %in% link_bins) |>
    filter(linked) |>
    mutate(cum.width = cumsum(width)) |>
    filter(cum.width <= size) |>
    pull(band)

  # Return bins from selected bands plus chr12q region
  lpsA_chromInfo |>
    filter(
      band %in% band_rois |
      (arm == "12q" & bin.start > 55e6 & bin.end < 105e6)
    ) |>
    pull(bin.id)
}


#' Summarize link types from get_links output
#'
#' Creates a summary table of structural variant link types with
#' metadata about the sample and filtering parameters.
#'
#' @param links_result Output from get_links() function
#' @param sample_id Sample identifier
#' @param tcn_filter Total copy number filter threshold used
#' @param pe_filter Paired-end read filter threshold used
#' @param subset_label Descriptive label for this subset (e.g., "Full", "EDC")
#' @return Data frame with link type counts as columns plus metadata
#'
#' @examples
#' links <- get_links(sample_data, 10, 3)
#' summarize_links(links, "WD0539P", 10, 3, "Full")
summarize_links <- function(links_result, sample_id, tcn_filter,
                           pe_filter, subset_label) {
  links_result$links |>
    count(Type) |>
    mutate(
      SampleID = sample_id,
      TCN.FILT = tcn_filter,
      PE.FILT = pe_filter,
      Subset = subset_label
    )
}