#' Summarise a nanopore run for QC handoff
#'
#' `qc_run_summary()` converts a normalised sequencing-summary table into a
#' stable one-row QC contract suitable for report cards, Grammateus reports,
#' and future Synoptikon handoff.
#'
#' @param x A `SequencingSummary` object, a sequencing-summary file path, or a
#'   data frame/tibble with floundeR-normalised sequencing-summary columns.
#' @param source_id Optional run, sample, file, or source identifier. When `x`
#'   is a path and `source_id` is not supplied, the path is used.
#' @param schema_version Schema version label to attach to the returned table.
#'
#' @return A one-row tibble with schema version, source identity, read/yield
#'   totals, pass/fail counts, length and Q-score summaries, channel count,
#'   observed time bounds, and barcode counts.
#'
#' @examples
#' summary_file <- flnDr("sequencing_summary.txt.bz2")
#' qc_run_summary(summary_file)
#'
#' @export
qc_run_summary <- function(
    x,
    source_id = NULL,
    schema_version = "flounder.qc_run_summary.v1") {
  seqsum <- .qc_run_summary_input(x)
  if (is.null(source_id) && is.character(x) && length(x) == 1) {
    source_id <- x
  }
  if (is.null(source_id)) {
    source_id <- NA_character_
  }

  required <- c(
    "channel", "start_time", "duration", "passes_filtering",
    "sequence_length_template", "mean_qscore_template")
  missing_required <- setdiff(required, names(seqsum))
  if (length(missing_required) > 0) {
    stop(
      "qc_run_summary requires sequencing-summary column(s): ",
      paste(missing_required, collapse = ", "),
      call. = FALSE)
  }
  if (!"barcode_arrangement" %in% names(seqsum)) {
    seqsum$barcode_arrangement <- NA_character_
  }

  read_count <- nrow(seqsum)
  length_values <- seqsum$sequence_length_template
  qscore_values <- seqsum$mean_qscore_template
  pass_values <- seqsum$passes_filtering
  passed_known <- !is.na(pass_values)

  passed_read_count <- sum(pass_values %in% TRUE, na.rm = TRUE)
  failed_read_count <- sum(pass_values %in% FALSE, na.rm = TRUE)
  pass_fraction <- if (any(passed_known)) {
    passed_read_count / sum(passed_known)
  } else {
    NA_real_
  }

  passed_bases <- .qc_sum(length_values[pass_values %in% TRUE])
  failed_bases <- .qc_sum(length_values[pass_values %in% FALSE])

  read_end <- seqsum$start_time + seqsum$duration
  first_read_start_time <- .qc_min(seqsum$start_time)
  last_read_start_time <- .qc_max(seqsum$start_time)
  run_duration_seconds <- .qc_max(read_end) - first_read_start_time
  if (is.na(first_read_start_time)) {
    run_duration_seconds <- NA_real_
  }

  barcode_values <- seqsum$barcode_arrangement
  classified_barcodes <- unique(
    barcode_values[
      !is.na(barcode_values) &
        barcode_values != "" &
        barcode_values != "unclassified"])

  tibble::tibble(
    schema_version = schema_version,
    source_id = source_id,
    read_count = read_count,
    passed_read_count = passed_read_count,
    failed_read_count = failed_read_count,
    pass_fraction = pass_fraction,
    total_bases = .qc_sum(length_values),
    passed_bases = passed_bases,
    failed_bases = failed_bases,
    mean_read_length = .qc_mean(length_values),
    median_read_length = .qc_median(length_values),
    n50_read_length = .qc_n50(length_values),
    mean_qscore = .qc_mean(qscore_values),
    median_qscore = .qc_median(qscore_values),
    channel_count = length(unique(stats::na.omit(seqsum$channel))),
    first_read_start_time = first_read_start_time,
    last_read_start_time = last_read_start_time,
    run_duration_seconds = run_duration_seconds,
    barcode_count = length(classified_barcodes),
    unclassified_read_count = sum(barcode_values == "unclassified", na.rm = TRUE)
  )
}

.qc_run_summary_input <- function(x) {
  if (inherits(x, "SequencingSummary")) {
    return(x$as_tibble())
  }

  if (is.character(x) && length(x) == 1) {
    return(SequencingSummary$new(x)$as_tibble())
  }

  if (is.data.frame(x)) {
    return(tibble::as_tibble(x))
  }

  stop(
    "qc_run_summary requires a SequencingSummary object, data frame, tibble, ",
    "or sequencing-summary file path.",
    call. = FALSE)
}

.qc_sum <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }

  sum(x, na.rm = TRUE)
}

.qc_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}

.qc_median <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }

  stats::median(x, na.rm = TRUE)
}

.qc_min <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }

  min(x, na.rm = TRUE)
}

.qc_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }

  max(x, na.rm = TRUE)
}

.qc_n50 <- function(x) {
  x <- stats::na.omit(x)
  if (length(x) == 0 || sum(x) <= 0) {
    return(NA_real_)
  }

  x <- sort(x, decreasing = TRUE)
  x[which(cumsum(x) >= sum(x) * 0.5)[[1]]]
}
