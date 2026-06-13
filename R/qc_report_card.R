#' Default report-card thresholds
#'
#' @return A named list of warn/fail thresholds used by `qc_report_card()`.
#'
#' @export
qc_report_card_thresholds <- function() {
  list(
    pass_fraction_min = c(warn = 0.85, fail = 0.70),
    mean_qscore_min = c(warn = 10, fail = 8),
    n50_read_length_min = c(warn = 1000, fail = 500),
    total_bases_min = c(warn = 1e6, fail = 1e5),
    channel_count_min = c(warn = 64, fail = 16),
    unclassified_fraction_max = c(warn = 0.10, fail = 0.25),
    barcode_max_fraction = c(warn = 0.60, fail = 0.80)
  )
}

#' Build a run-level QC report card
#'
#' `qc_report_card()` evaluates a run summary against pass/warn/fail threshold
#' checks. It is intentionally limited to sequencing-summary-derived evidence
#' until POD5, BAM/Bamana, and Porkchop contracts are available.
#'
#' @param x A `SequencingSummary` object, a sequencing-summary file path, a
#'   normalised sequencing-summary data frame/tibble, or a `qc_run_summary()`
#'   one-row tibble.
#' @param thresholds Named threshold list. Defaults to
#'   `qc_report_card_thresholds()`.
#' @param barcode_composition Optional `qc_barcode_composition()` tibble. When
#'   omitted and `x` is sequencing-summary-like, it is derived automatically.
#' @param schema_version Schema version label to attach to the returned table.
#'
#' @return A tibble with one row per report-card check and columns:
#'   `schema_version`, `check_id`, `check_label`, `status`, `observed_value`,
#'   `warn_threshold`, `fail_threshold`, `comparator`, and `details`.
#'
#' @examples
#' summary_file <- flnDr("sequencing_summary.txt.bz2")
#' qc_report_card(summary_file)
#'
#' @export
qc_report_card <- function(
    x,
    thresholds = qc_report_card_thresholds(),
    barcode_composition = NULL,
    schema_version = "flounder.qc_report_card.v1") {
  run_summary <- .qc_report_card_summary(x)

  if (is.null(barcode_composition) && !.qc_is_run_summary(x)) {
    barcode_composition <- tryCatch(
      qc_barcode_composition(x),
      error = function(e) NULL)
  }

  unclassified_fraction <- .qc_unclassified_fraction(run_summary)
  barcode_max_fraction <- .qc_barcode_max_fraction(barcode_composition)

  .qc_bind_rows(list(
    .qc_threshold_check(
      schema_version,
      check_id = "pass_fraction",
      check_label = "Pass read fraction",
      observed_value = run_summary$pass_fraction[[1]],
      threshold = thresholds$pass_fraction_min,
      comparator = "minimum"),
    .qc_threshold_check(
      schema_version,
      check_id = "mean_qscore",
      check_label = "Mean read Q-score",
      observed_value = run_summary$mean_qscore[[1]],
      threshold = thresholds$mean_qscore_min,
      comparator = "minimum"),
    .qc_threshold_check(
      schema_version,
      check_id = "n50_read_length",
      check_label = "Read length N50",
      observed_value = run_summary$n50_read_length[[1]],
      threshold = thresholds$n50_read_length_min,
      comparator = "minimum"),
    .qc_threshold_check(
      schema_version,
      check_id = "total_bases",
      check_label = "Total base yield",
      observed_value = run_summary$total_bases[[1]],
      threshold = thresholds$total_bases_min,
      comparator = "minimum"),
    .qc_threshold_check(
      schema_version,
      check_id = "channel_count",
      check_label = "Active channel count",
      observed_value = run_summary$channel_count[[1]],
      threshold = thresholds$channel_count_min,
      comparator = "minimum"),
    .qc_threshold_check(
      schema_version,
      check_id = "unclassified_fraction",
      check_label = "Unclassified barcode fraction",
      observed_value = unclassified_fraction,
      threshold = thresholds$unclassified_fraction_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "barcode_max_fraction",
      check_label = "Largest barcode read fraction",
      observed_value = barcode_max_fraction,
      threshold = thresholds$barcode_max_fraction,
      comparator = "maximum")
  ))
}

.qc_report_card_summary <- function(x) {
  if (.qc_is_run_summary(x)) {
    return(tibble::as_tibble(x))
  }

  qc_run_summary(x)
}

.qc_is_run_summary <- function(x) {
  is.data.frame(x) &&
    all(c(
      "schema_version", "read_count", "pass_fraction", "total_bases",
      "mean_qscore", "n50_read_length", "channel_count",
      "unclassified_read_count") %in% names(x)) &&
    length(x$schema_version) >= 1 &&
    grepl("^flounder[.]qc_run_summary[.]v", x$schema_version[[1]])
}

.qc_threshold_check <- function(
    schema_version,
    check_id,
    check_label,
    observed_value,
    threshold,
    comparator) {
  threshold <- .qc_threshold_pair(threshold, check_id)
  status <- .qc_threshold_status(observed_value, threshold, comparator)

  tibble::tibble(
    schema_version = schema_version,
    check_id = check_id,
    check_label = check_label,
    status = status,
    observed_value = as.numeric(observed_value),
    warn_threshold = as.numeric(threshold[["warn"]]),
    fail_threshold = as.numeric(threshold[["fail"]]),
    comparator = comparator,
    details = .qc_threshold_details(
      observed_value, threshold, comparator, status)
  )
}

.qc_threshold_pair <- function(threshold, check_id) {
  if (is.null(threshold) ||
      !all(c("warn", "fail") %in% names(threshold)) ||
      any(is.na(threshold[c("warn", "fail")]))) {
    stop(
      "Threshold for `", check_id, "` must contain named warn and fail values.",
      call. = FALSE)
  }

  threshold[c("warn", "fail")]
}

.qc_threshold_status <- function(observed_value, threshold, comparator) {
  if (length(observed_value) == 0 || is.na(observed_value)) {
    return("warn")
  }

  if (comparator == "minimum") {
    if (observed_value < threshold[["fail"]]) {
      return("fail")
    }
    if (observed_value < threshold[["warn"]]) {
      return("warn")
    }
    return("pass")
  }

  if (comparator == "maximum") {
    if (observed_value > threshold[["fail"]]) {
      return("fail")
    }
    if (observed_value > threshold[["warn"]]) {
      return("warn")
    }
    return("pass")
  }

  stop("Unsupported threshold comparator: ", comparator, call. = FALSE)
}

.qc_threshold_details <- function(
    observed_value,
    threshold,
    comparator,
    status) {
  if (length(observed_value) == 0 || is.na(observed_value)) {
    return("Observed value is missing; review source data availability.")
  }

  if (status == "pass") {
    return("Observed value is within the pass range.")
  }

  if (status == "warn") {
    return(paste0(
      "Observed value crossed the warning threshold (",
      comparator, " warn=", threshold[["warn"]], ")."))
  }

  paste0(
    "Observed value crossed the failure threshold (",
    comparator, " fail=", threshold[["fail"]], ").")
}

.qc_unclassified_fraction <- function(run_summary) {
  read_count <- run_summary$read_count[[1]]
  if (is.na(read_count) || read_count == 0) {
    return(NA_real_)
  }

  run_summary$unclassified_read_count[[1]] / read_count
}

.qc_barcode_max_fraction <- function(barcode_composition) {
  if (is.null(barcode_composition) ||
      !"barcode_arrangement" %in% names(barcode_composition) ||
      !"read_fraction" %in% names(barcode_composition)) {
    return(NA_real_)
  }

  classified <- barcode_composition[
    !is.na(barcode_composition$barcode_arrangement) &
      barcode_composition$barcode_arrangement != "unclassified",
    , drop = FALSE]
  if (nrow(classified) == 0 || all(is.na(classified$read_fraction))) {
    return(NA_real_)
  }

  max(classified$read_fraction, na.rm = TRUE)
}
