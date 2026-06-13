#' Build tidy QC data from sequencing-summary inputs
#'
#' These helpers return stable, schema-versioned tibbles for the first
#' floundeR QC plot and reporting families. They accept the same input shapes
#' as `qc_run_summary()`: a `SequencingSummary` object, a sequencing-summary
#' file path, or a normalised data frame/tibble.
#'
#' @param x A `SequencingSummary` object, a sequencing-summary file path, or a
#'   data frame/tibble with floundeR-normalised sequencing-summary columns.
#' @param resolution_minutes Temporal bin size in minutes.
#' @param bins Number of distribution bins to produce.
#'
#' @return A tibble whose columns are documented by each function.
NULL

#' Yield Over Time QC Contract
#'
#' @inheritParams qc_run_summary
#' @param resolution_minutes Temporal bin size in minutes.
#'
#' @return A tibble with schema version, bin bounds, pass/fail state, read
#'   count, base yield, and cumulative read/base yield.
#'
#' @export
qc_yield_over_time <- function(
    x,
    resolution_minutes = 15,
    schema_version = "flounder.qc_yield_over_time.v1") {
  if (!is.numeric(resolution_minutes) || length(resolution_minutes) != 1 ||
      is.na(resolution_minutes) || resolution_minutes <= 0) {
    stop("resolution_minutes must be a positive numeric scalar.", call. = FALSE)
  }

  seqsum <- .qc_prepare_seqsum(
    x,
    c(
      "start_time", "passes_filtering", "sequence_length_template"))

  resolution_seconds <- resolution_minutes * 60
  bin_start_seconds <- floor(seqsum$start_time / resolution_seconds) *
    resolution_seconds
  seqsum$.bin_start_seconds <- bin_start_seconds
  seqsum$.bin_end_seconds <- bin_start_seconds + resolution_seconds

  grouped <- .qc_group_counts(
    seqsum,
    c(".bin_start_seconds", ".bin_end_seconds", "passes_filtering"))
  grouped <- grouped[order(
    grouped$.bin_start_seconds,
    is.na(grouped$passes_filtering),
    grouped$passes_filtering), , drop = FALSE]

  grouped$cumulative_read_count <- stats::ave(
    grouped$read_count,
    grouped$passes_filtering,
    FUN = cumsum)
  grouped$cumulative_bases <- stats::ave(
    grouped$bases,
    grouped$passes_filtering,
    FUN = cumsum)

  tibble::tibble(
    schema_version = schema_version,
    bin_start_seconds = grouped$.bin_start_seconds,
    bin_end_seconds = grouped$.bin_end_seconds,
    bin_start_minutes = grouped$.bin_start_seconds / 60,
    bin_end_minutes = grouped$.bin_end_seconds / 60,
    passes_filtering = grouped$passes_filtering,
    read_count = grouped$read_count,
    bases = grouped$bases,
    cumulative_read_count = grouped$cumulative_read_count,
    cumulative_bases = grouped$cumulative_bases
  )
}

#' Read Length Distribution QC Contract
#'
#' @inheritParams qc_run_summary
#' @param bins Number of read-length bins to produce.
#'
#' @return A tibble with schema version, read-length bin bounds, pass/fail
#'   state, read count, and base yield.
#'
#' @export
qc_read_length_distribution <- function(
    x,
    bins = 20,
    schema_version = "flounder.qc_read_length_distribution.v1") {
  .qc_distribution(
    x = x,
    column = "sequence_length_template",
    bins = bins,
    schema_version = schema_version,
    prefix = "read_length")
}

#' Quality Distribution QC Contract
#'
#' @inheritParams qc_run_summary
#' @param bins Number of Q-score bins to produce.
#'
#' @return A tibble with schema version, Q-score bin bounds, pass/fail state,
#'   read count, and base yield.
#'
#' @export
qc_quality_distribution <- function(
    x,
    bins = 20,
    schema_version = "flounder.qc_quality_distribution.v1") {
  .qc_distribution(
    x = x,
    column = "mean_qscore_template",
    bins = bins,
    schema_version = schema_version,
    prefix = "qscore")
}

#' Channel Density QC Contract
#'
#' @inheritParams qc_run_summary
#'
#' @return A tibble with schema version, channel, read/base yield, pass/fail
#'   read counts, and pass fraction.
#'
#' @export
qc_channel_density <- function(
    x,
    schema_version = "flounder.qc_channel_density.v1") {
  seqsum <- .qc_prepare_seqsum(
    x,
    c("channel", "passes_filtering", "sequence_length_template"))
  split_rows <- split(seqsum, seqsum$channel, drop = TRUE)
  rows <- lapply(split_rows, function(group) {
    passed <- group$passes_filtering %in% TRUE
    failed <- group$passes_filtering %in% FALSE
    passed_known <- !is.na(group$passes_filtering)
    passed_read_count <- sum(passed)
    tibble::tibble(
      schema_version = schema_version,
      channel = group$channel[[1]],
      read_count = nrow(group),
      bases = .qc_sum(group$sequence_length_template),
      passed_read_count = passed_read_count,
      failed_read_count = sum(failed),
      pass_fraction = if (any(passed_known)) {
        passed_read_count / sum(passed_known)
      } else {
        NA_real_
      }
    )
  })

  .qc_bind_rows(rows)
}

#' Barcode Composition QC Contract
#'
#' @inheritParams qc_run_summary
#'
#' @return A tibble with schema version, barcode arrangement, read/base yield,
#'   pass/fail read counts, and read/base fractions.
#'
#' @export
qc_barcode_composition <- function(
    x,
    schema_version = "flounder.qc_barcode_composition.v1") {
  seqsum <- .qc_prepare_seqsum(
    x,
    c("passes_filtering", "sequence_length_template"),
    optional = "barcode_arrangement")
  if (!"barcode_arrangement" %in% names(seqsum)) {
    seqsum$barcode_arrangement <- "unclassified"
  }
  seqsum$barcode_arrangement[is.na(seqsum$barcode_arrangement) |
    seqsum$barcode_arrangement == ""] <- "unclassified"

  total_reads <- nrow(seqsum)
  total_bases <- .qc_sum(seqsum$sequence_length_template)
  split_rows <- split(seqsum, seqsum$barcode_arrangement, drop = TRUE)
  rows <- lapply(split_rows, function(group) {
    bases <- .qc_sum(group$sequence_length_template)
    tibble::tibble(
      schema_version = schema_version,
      barcode_arrangement = group$barcode_arrangement[[1]],
      read_count = nrow(group),
      bases = bases,
      passed_read_count = sum(group$passes_filtering %in% TRUE),
      failed_read_count = sum(group$passes_filtering %in% FALSE),
      read_fraction = if (total_reads > 0) {
        nrow(group) / total_reads
      } else {
        NA_real_
      },
      bases_fraction = if (!is.na(total_bases) && total_bases > 0) {
        bases / total_bases
      } else {
        NA_real_
      }
    )
  })

  .qc_bind_rows(rows)
}

.qc_prepare_seqsum <- function(x, required, optional = character()) {
  seqsum <- .qc_run_summary_input(x)
  missing_required <- setdiff(required, names(seqsum))
  if (length(missing_required) > 0) {
    stop(
      "QC contract requires sequencing-summary column(s): ",
      paste(missing_required, collapse = ", "),
      call. = FALSE)
  }

  seqsum[, unique(c(required, optional))[unique(c(required, optional)) %in%
    names(seqsum)], drop = FALSE]
}

.qc_distribution <- function(x, column, bins, schema_version, prefix) {
  if (!is.numeric(bins) || length(bins) != 1 || is.na(bins) || bins < 1) {
    stop("bins must be a positive numeric scalar.", call. = FALSE)
  }
  bins <- as.integer(bins)
  seqsum <- .qc_prepare_seqsum(
    x,
    c(column, "passes_filtering", "sequence_length_template"))
  values <- seqsum[[column]]
  if (all(is.na(values))) {
    stop("Distribution column contains only missing values.", call. = FALSE)
  }

  max_value <- max(values, na.rm = TRUE)
  breaks <- seq(0, max_value, length.out = bins + 1)
  if (length(unique(breaks)) == 1) {
    breaks <- c(0, max_value + 1)
  }

  seqsum$.bin <- cut(
    values,
    breaks = breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = FALSE)
  seqsum$.bin_start <- breaks[seqsum$.bin]
  seqsum$.bin_end <- breaks[seqsum$.bin + 1]

  grouped <- .qc_group_counts(
    seqsum,
    c(".bin_start", ".bin_end", "passes_filtering"))
  grouped <- grouped[order(
    grouped$.bin_start,
    is.na(grouped$passes_filtering),
    grouped$passes_filtering), , drop = FALSE]

  result <- tibble::tibble(
    schema_version = schema_version,
    passes_filtering = grouped$passes_filtering,
    read_count = grouped$read_count,
    bases = grouped$bases
  )
  result[[paste0(prefix, "_bin_start")]] <- grouped$.bin_start
  result[[paste0(prefix, "_bin_end")]] <- grouped$.bin_end
  result[, c(
    "schema_version",
    paste0(prefix, "_bin_start"),
    paste0(prefix, "_bin_end"),
    "passes_filtering",
    "read_count",
    "bases")]
}

.qc_group_counts <- function(seqsum, group_cols) {
  rows <- split(
    seqsum,
    interaction(seqsum[, group_cols, drop = FALSE], drop = TRUE),
    drop = TRUE)
  result <- lapply(rows, function(group) {
    keys <- group[1, group_cols, drop = FALSE]
    keys$read_count <- nrow(group)
    keys$bases <- .qc_sum(group$sequence_length_template)
    keys
  })

  .qc_bind_rows(result)
}

.qc_bind_rows <- function(rows) {
  rows <- rows[lengths(rows) > 0]
  if (length(rows) == 0) {
    return(tibble::tibble())
  }

  tibble::as_tibble(do.call(rbind, rows))
}
