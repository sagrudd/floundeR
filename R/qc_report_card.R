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

#' Default BAM report-card thresholds
#'
#' @return A named list of warn/fail thresholds used by
#'   `bam_qc_report_card()`.
#'
#' @export
bam_qc_report_card_thresholds <- function() {
  list(
    mapping_fraction_min = c(warn = 0.80, fail = 0.50),
    duplicate_fraction_max = c(warn = 0.20, fail = 0.40),
    qc_fail_fraction_max = c(warn = 0.05, fail = 0.15),
    mapq_zero_fraction_max = c(warn = 0.10, fail = 0.25),
    missing_or_unusable_index_max = c(warn = 0, fail = 0),
    stale_index_max = c(warn = 0, fail = 0),
    sorting_mismatch_max = c(warn = 0, fail = 0),
    missing_expected_tags_max = c(warn = 0, fail = 0),
    eof_absence_max = c(warn = 0, fail = 0),
    validation_finding_count_max = c(warn = 0, fail = 0),
    provenance_anomaly_count_max = c(warn = 0, fail = 0)
  )
}

#' Default library-preparation report-card thresholds
#'
#' @return A named list of warn/fail thresholds used by
#'   `library_preparation_report_card()`.
#'
#' @export
library_preparation_report_card_thresholds <- function() {
  list(
    unexpected_kit_max = c(warn = 0, fail = 0),
    adapter_burden_max = c(warn = 0.20, fail = 0.50),
    barcode_ambiguity_fraction_max = c(warn = 0.05, fail = 0.20),
    cdna_incomplete_fraction_max = c(warn = 0.20, fail = 0.50),
    kit_support_risk_max = c(warn = 0, fail = 1)
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

#' Build a BAM QC report card
#'
#' `bam_qc_report_card()` evaluates Bamana-derived BAM evidence against
#' pass/warn/fail checks for alignment QC and reporting. It consumes the
#' R-native objects returned by `bam_summary()`, `bam_check_index()`,
#' `bam_check_sort()`, `bam_check_tag()`, `bam_check_eof()`, and
#' `bam_validate()`; it does not call Bamana directly.
#'
#' @param summary Optional `bam_summary()` result.
#' @param index Optional `bam_check_index()` result.
#' @param mapping Optional `bam_check_map()` result used as a fallback for
#'   mapping-fraction evidence when `summary` is absent.
#' @param sorting Optional `bam_check_sort()` result.
#' @param tags Optional `bam_check_tag()` data frame or list of such data
#'   frames.
#' @param eof Optional `bam_check_eof()` data frame.
#' @param validation Optional `bam_validate()` result.
#' @param provenance_anomalies Optional numeric anomaly count or data frame of
#'   provenance findings from a future Bamana provenance/forensic surface.
#' @param expected_tags Character vector of BAM aux tags expected to be present.
#' @param thresholds Named threshold list. Defaults to
#'   `bam_qc_report_card_thresholds()`.
#' @param schema_version Schema version label to attach to the returned table.
#'
#' @return A tibble using the standard report-card schema:
#'   `schema_version`, `check_id`, `check_label`, `status`, `observed_value`,
#'   `warn_threshold`, `fail_threshold`, `comparator`, and `details`.
#'
#' @examples
#' \dontrun{
#' bam <- bam_summary("reads.bam")
#' bam_qc_report_card(summary = bam)
#' }
#'
#' @export
bam_qc_report_card <- function(
    summary = NULL,
    index = NULL,
    mapping = NULL,
    sorting = NULL,
    tags = NULL,
    eof = NULL,
    validation = NULL,
    provenance_anomalies = NULL,
    expected_tags = character(),
    thresholds = bam_qc_report_card_thresholds(),
    schema_version = "flounder.bam_qc_report_card.v1") {
  if (!is.character(expected_tags) || any(is.na(expected_tags))) {
    stop("`expected_tags` must be a character vector without NA values.",
         call. = FALSE)
  }

  .qc_bind_rows(list(
    .qc_threshold_check(
      schema_version,
      check_id = "bam_mapping_fraction",
      check_label = "BAM mapped-read fraction",
      observed_value = .bam_card_mapping_fraction(summary, mapping),
      threshold = thresholds$mapping_fraction_min,
      comparator = "minimum"),
    .qc_threshold_check(
      schema_version,
      check_id = "bam_duplicate_fraction",
      check_label = "BAM duplicate-record fraction",
      observed_value = .bam_card_summary_fraction(summary, "fraction_duplicate"),
      threshold = thresholds$duplicate_fraction_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "bam_qc_fail_fraction",
      check_label = "BAM QC-fail-record fraction",
      observed_value = .bam_card_summary_fraction(summary, "fraction_qc_fail"),
      threshold = thresholds$qc_fail_fraction_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "bam_mapq_zero_fraction",
      check_label = "BAM MAPQ-zero fraction",
      observed_value = .bam_card_mapq_zero_fraction(summary),
      threshold = thresholds$mapq_zero_fraction_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "bam_missing_or_unusable_index",
      check_label = "Missing or unusable BAM index",
      observed_value = .bam_card_missing_or_unusable_index(index),
      threshold = thresholds$missing_or_unusable_index_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "bam_stale_index",
      check_label = "Stale BAM index",
      observed_value = .bam_card_stale_index(index),
      threshold = thresholds$stale_index_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "bam_sorting_mismatch",
      check_label = "BAM sorting mismatch",
      observed_value = .bam_card_sorting_mismatch(sorting),
      threshold = thresholds$sorting_mismatch_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "bam_missing_expected_tags",
      check_label = "Missing expected BAM aux tags",
      observed_value = .bam_card_missing_expected_tags(tags, expected_tags),
      threshold = thresholds$missing_expected_tags_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "bam_eof_absence",
      check_label = "Missing BGZF EOF marker",
      observed_value = .bam_card_eof_absence(eof),
      threshold = thresholds$eof_absence_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "bam_validation_findings",
      check_label = "BAM validation findings",
      observed_value = .bam_card_validation_finding_count(validation),
      threshold = thresholds$validation_finding_count_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "bam_provenance_anomalies",
      check_label = "BAM provenance anomalies",
      observed_value = .bam_card_provenance_anomaly_count(provenance_anomalies),
      threshold = thresholds$provenance_anomaly_count_max,
      comparator = "maximum")
  ))
}

#' Build a library-preparation QC report card
#'
#' `library_preparation_report_card()` evaluates Porkchop-derived
#' library-preparation evidence against pass/warn/fail checks for reporting and
#' review. It consumes R-native evidence tables returned by
#' `library_kit_candidates()`, `library_adapter_primer_evidence()`,
#' `library_barcode_evidence()`, and `library_cdna_primer_evidence()`; it does
#' not call Porkchop directly.
#'
#' @param kit_candidates Optional `library_kit_candidates()` data frame.
#' @param adapter_primer Optional `library_adapter_primer_evidence()` data
#'   frame.
#' @param barcode Optional `library_barcode_evidence()` data frame.
#' @param cdna Optional `library_cdna_primer_evidence()` data frame.
#' @param expected_kit_id Optional scalar expected kit id from run metadata or
#'   a user-supplied review context.
#' @param thresholds Named threshold list. Defaults to
#'   `library_preparation_report_card_thresholds()`.
#' @param schema_version Schema version label to attach to the returned table.
#'
#' @return A tibble using the standard report-card schema:
#'   `schema_version`, `check_id`, `check_label`, `status`, `observed_value`,
#'   `warn_threshold`, `fail_threshold`, `comparator`, and `details`.
#'
#' @examples
#' kit_candidates <- data.frame(
#'   kit_id = "SQK-LSK114",
#'   normalized_score = 0.9,
#'   support_level = "supported",
#'   lifecycle_status = "current"
#' )
#' library_preparation_report_card(
#'   kit_candidates = kit_candidates,
#'   expected_kit_id = "SQK-LSK114"
#' )
#'
#' @export
library_preparation_report_card <- function(
    kit_candidates = NULL,
    adapter_primer = NULL,
    barcode = NULL,
    cdna = NULL,
    expected_kit_id = NULL,
    thresholds = library_preparation_report_card_thresholds(),
    schema_version = "flounder.library_preparation_report_card.v1") {
  .library_card_validate_optional_frame(kit_candidates, "kit_candidates")
  .library_card_validate_optional_frame(adapter_primer, "adapter_primer")
  .library_card_validate_optional_frame(barcode, "barcode")
  .library_card_validate_optional_frame(cdna, "cdna")
  expected_kit_id <- .library_card_expected_kit_id(expected_kit_id)

  top_candidate <- .library_card_top_candidate(kit_candidates)

  .qc_bind_rows(list(
    .qc_threshold_check(
      schema_version,
      check_id = "library_unexpected_kit",
      check_label = "Unexpected library kit",
      observed_value = .library_card_unexpected_kit(
        top_candidate, expected_kit_id),
      threshold = thresholds$unexpected_kit_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "library_adapter_burden",
      check_label = "Adapter motif burden",
      observed_value = .library_card_adapter_burden(
        adapter_primer, barcode, cdna),
      threshold = thresholds$adapter_burden_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "library_barcode_ambiguity",
      check_label = "Barcode ambiguity fraction",
      observed_value = .library_card_barcode_ambiguity(barcode),
      threshold = thresholds$barcode_ambiguity_fraction_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "library_cdna_incomplete_fraction",
      check_label = "cDNA partial or unclassified fraction",
      observed_value = .library_card_cdna_incomplete_fraction(cdna),
      threshold = thresholds$cdna_incomplete_fraction_max,
      comparator = "maximum"),
    .qc_threshold_check(
      schema_version,
      check_id = "library_kit_support_risk",
      check_label = "Unsupported or partial kit support level",
      observed_value = .library_card_support_level_risk(top_candidate),
      threshold = thresholds$kit_support_risk_max,
      comparator = "maximum")
  ))
}

.qc_report_card_summary <- function(x) {
  if (.qc_is_run_summary(x)) {
    return(tibble::as_tibble(x))
  }

  qc_run_summary(x)
}

.bam_card_summary_fraction <- function(summary, column) {
  if (is.null(summary) || is.null(summary$fractions) ||
      !column %in% names(summary$fractions)) {
    return(NA_real_)
  }

  fractions <- summary$fractions
  if ("scope" %in% names(fractions) && any(fractions$scope == "full_file")) {
    fractions <- fractions[fractions$scope == "full_file", , drop = FALSE]
  }

  as.numeric(fractions[[column]][[1]])
}

.bam_card_mapping_fraction <- function(summary, mapping) {
  fraction <- .bam_card_summary_fraction(summary, "fraction_mapped")
  if (!is.na(fraction)) {
    return(fraction)
  }

  if (is.null(mapping) || is.null(mapping$summary)) {
    return(NA_real_)
  }

  mapped <- .bam_card_first_numeric(mapping$summary, "total_mapped_reads")
  unmapped <- .bam_card_first_numeric(mapping$summary, "total_unmapped_reads")
  if (is.na(mapped) || is.na(unmapped) || mapped + unmapped == 0) {
    mapped <- .bam_card_first_numeric(mapping$summary, "mapped_records_observed")
    unmapped <- .bam_card_first_numeric(mapping$summary, "unmapped_records_observed")
  }
  if (is.na(mapped) || is.na(unmapped) || mapped + unmapped == 0) {
    return(NA_real_)
  }

  mapped / (mapped + unmapped)
}

.bam_card_mapq_zero_fraction <- function(summary) {
  if (is.null(summary) || is.null(summary$mapq) || is.null(summary$counts)) {
    return(NA_real_)
  }

  zero_count <- .bam_card_first_numeric(summary$mapq, "zero_count")
  records <- .bam_card_first_numeric(summary$counts, "records_examined")
  if (is.na(zero_count) || is.na(records) || records == 0) {
    return(NA_real_)
  }

  zero_count / records
}

.bam_card_missing_or_unusable_index <- function(index) {
  if (is.null(index) || is.null(index$index)) {
    return(NA_real_)
  }

  present <- .bam_card_first_logical(index$index, "present")
  usable <- .bam_card_first_logical(index$index, "usable")
  if (is.na(present) || is.na(usable)) {
    return(NA_real_)
  }

  as.numeric(!present || !usable)
}

.bam_card_stale_index <- function(index) {
  if (is.null(index) || is.null(index$index)) {
    return(NA_real_)
  }

  stale <- .bam_card_first_bool_label(index$index, "stale")
  bam_newer <- .bam_card_first_bool_label(index$index, "bam_newer_than_index")
  compatibility <- .bam_card_first_character(index$index, "compatibility")
  if (is.na(stale) && is.na(bam_newer) && is.na(compatibility)) {
    return(NA_real_)
  }

  as.numeric(
    identical(stale, TRUE) ||
      identical(bam_newer, TRUE) ||
      identical(compatibility, "stale"))
}

.bam_card_sorting_mismatch <- function(sorting) {
  if (is.null(sorting)) {
    return(NA_real_)
  }

  header_match <- .bam_card_first_bool_label(sorting, "header_matches_observation")
  appears_sorted <- .bam_card_first_bool_label(sorting, "appears_sorted")
  violation <- .bam_card_first_character(sorting, "first_violation_reason")
  if (is.na(header_match) && is.na(appears_sorted) && is.na(violation)) {
    return(NA_real_)
  }

  as.numeric(
    identical(header_match, FALSE) ||
      identical(appears_sorted, FALSE) ||
      !is.na(violation))
}

.bam_card_missing_expected_tags <- function(tags, expected_tags) {
  expected_tags <- unique(expected_tags)
  if (length(expected_tags) == 0) {
    return(0)
  }
  if (is.null(tags)) {
    return(NA_real_)
  }

  tag_table <- .bam_card_tag_table(tags)
  if (nrow(tag_table) == 0 || !"tag" %in% names(tag_table) ||
      !"tag_found" %in% names(tag_table)) {
    return(NA_real_)
  }

  found <- unique(tag_table$tag[tag_table$tag %in% expected_tags &
    tag_table$tag_found %in% TRUE])
  length(setdiff(expected_tags, found))
}

.bam_card_eof_absence <- function(eof) {
  if (is.null(eof) || !"complete" %in% names(eof)) {
    return(NA_real_)
  }

  complete <- eof$complete[[1]]
  if (is.na(complete)) {
    return(NA_real_)
  }

  as.numeric(!complete)
}

.bam_card_validation_finding_count <- function(validation) {
  if (is.null(validation)) {
    return(NA_real_)
  }
  if (!is.null(validation$findings)) {
    return(nrow(validation$findings))
  }
  if (is.null(validation$summary)) {
    return(NA_real_)
  }

  sum(
    .bam_card_first_numeric(validation$summary, "errors"),
    .bam_card_first_numeric(validation$summary, "warnings"),
    .bam_card_first_numeric(validation$summary, "infos"),
    na.rm = TRUE)
}

.bam_card_provenance_anomaly_count <- function(provenance_anomalies) {
  if (is.null(provenance_anomalies)) {
    return(NA_real_)
  }
  if (is.numeric(provenance_anomalies) && length(provenance_anomalies) == 1L) {
    return(as.numeric(provenance_anomalies))
  }
  if (is.data.frame(provenance_anomalies)) {
    return(nrow(provenance_anomalies))
  }

  NA_real_
}

.bam_card_tag_table <- function(tags) {
  if (is.data.frame(tags)) {
    return(tags)
  }
  if (is.list(tags)) {
    frames <- tags[vapply(tags, is.data.frame, logical(1))]
    if (length(frames) == 0) {
      return(data.frame())
    }
    return(tibble::as_tibble(do.call(rbind, frames)))
  }

  data.frame()
}

.bam_card_first_numeric <- function(data, column) {
  if (is.null(data) || !column %in% names(data) || length(data[[column]]) == 0) {
    return(NA_real_)
  }

  as.numeric(data[[column]][[1]])
}

.bam_card_first_character <- function(data, column) {
  if (is.null(data) || !column %in% names(data) || length(data[[column]]) == 0) {
    return(NA_character_)
  }

  value <- as.character(data[[column]][[1]])
  if (is.na(value) || identical(value, "")) {
    return(NA_character_)
  }
  value
}

.bam_card_first_logical <- function(data, column) {
  if (is.null(data) || !column %in% names(data) || length(data[[column]]) == 0) {
    return(NA)
  }

  value <- data[[column]][[1]]
  if (is.na(value)) {
    return(NA)
  }
  as.logical(value)
}

.bam_card_first_bool_label <- function(data, column) {
  value <- .bam_card_first_character(data, column)
  if (is.na(value)) {
    return(NA)
  }
  if (tolower(value) == "true") {
    return(TRUE)
  }
  if (tolower(value) == "false") {
    return(FALSE)
  }

  NA
}

.library_card_validate_optional_frame <- function(data, name) {
  if (!is.null(data) && !is.data.frame(data)) {
    stop("`", name, "` must be a data frame or NULL.", call. = FALSE)
  }
}

.library_card_expected_kit_id <- function(expected_kit_id) {
  if (is.null(expected_kit_id)) {
    return(NULL)
  }
  if (!is.character(expected_kit_id) || length(expected_kit_id) != 1L ||
      is.na(expected_kit_id) || identical(expected_kit_id, "")) {
    stop("`expected_kit_id` must be a non-empty character scalar.", call. = FALSE)
  }

  expected_kit_id
}

.library_card_top_candidate <- function(kit_candidates) {
  if (is.null(kit_candidates) || nrow(kit_candidates) == 0 ||
      !"kit_id" %in% names(kit_candidates)) {
    return(NULL)
  }

  candidates <- tibble::as_tibble(kit_candidates)
  order_score <- rep(0, nrow(candidates))
  if ("normalized_score" %in% names(candidates)) {
    order_score <- suppressWarnings(as.numeric(candidates$normalized_score))
  } else if ("score" %in% names(candidates)) {
    order_score <- suppressWarnings(as.numeric(candidates$score))
  }
  order_score[is.na(order_score)] <- -Inf

  candidates[order(order_score, decreasing = TRUE), , drop = FALSE][1, ]
}

.library_card_unexpected_kit <- function(top_candidate, expected_kit_id) {
  if (is.null(expected_kit_id)) {
    return(0)
  }
  if (is.null(top_candidate) || !"kit_id" %in% names(top_candidate)) {
    return(NA_real_)
  }

  observed_kit <- as.character(top_candidate$kit_id[[1]])
  if (is.na(observed_kit) || identical(observed_kit, "")) {
    return(NA_real_)
  }

  as.numeric(!identical(observed_kit, expected_kit_id))
}

.library_card_adapter_burden <- function(adapter_primer, barcode, cdna) {
  if (is.null(adapter_primer) || nrow(adapter_primer) == 0 ||
      !"read_id" %in% names(adapter_primer)) {
    return(NA_real_)
  }

  denominator <- .library_card_read_count(adapter_primer, barcode, cdna)
  if (is.na(denominator) || denominator == 0) {
    return(NA_real_)
  }

  adapter_reads <- unique(adapter_primer$read_id[
    .library_card_row_matches(adapter_primer, "adapter")
  ])
  length(stats::na.omit(adapter_reads)) / denominator
}

.library_card_barcode_ambiguity <- function(barcode) {
  if (is.null(barcode) || nrow(barcode) == 0 || !"read_id" %in% names(barcode)) {
    return(NA_real_)
  }

  read_ids <- unique(stats::na.omit(barcode$read_id))
  if (length(read_ids) == 0) {
    return(NA_real_)
  }

  ambiguous <- vapply(read_ids, function(read_id) {
    rows <- barcode[barcode$read_id == read_id, , drop = FALSE]
    barcode_rows <- .library_card_row_matches(rows, "barcode")
    flank_rows <- .library_card_row_matches(rows, "flank")
    barcode_names <- unique(stats::na.omit(.library_card_motif_names(
      rows[barcode_rows, , drop = FALSE])))

    length(barcode_names) > 1L || (length(barcode_names) == 0L && any(flank_rows))
  }, logical(1))

  mean(ambiguous)
}

.library_card_cdna_incomplete_fraction <- function(cdna) {
  if (is.null(cdna) || nrow(cdna) == 0) {
    return(NA_real_)
  }

  classified <- rep(NA, nrow(cdna))
  if ("classified" %in% names(cdna)) {
    classified <- .library_card_as_logical(cdna$classified)
  }

  full_length <- rep(NA, nrow(cdna))
  if ("full_length" %in% names(cdna)) {
    full_length <- .library_card_as_logical(cdna$full_length)
  }

  class_label <- rep(NA_character_, nrow(cdna))
  if ("class" %in% names(cdna)) {
    class_label <- tolower(as.character(cdna$class))
  }

  incomplete <- classified %in% FALSE |
    full_length %in% FALSE |
    (!is.na(class_label) & class_label != "full_length")

  if (all(is.na(classified)) && all(is.na(full_length)) &&
      all(is.na(class_label))) {
    return(NA_real_)
  }

  mean(incomplete, na.rm = TRUE)
}

.library_card_support_level_risk <- function(top_candidate) {
  if (is.null(top_candidate)) {
    return(NA_real_)
  }

  support_level <- .bam_card_first_character(top_candidate, "support_level")
  lifecycle_status <- .bam_card_first_character(top_candidate, "lifecycle_status")
  validation_status <- .bam_card_first_character(top_candidate, "validation_status")

  risk <- 0
  labels <- tolower(stats::na.omit(c(
    support_level, lifecycle_status, validation_status)))
  if (length(labels) == 0) {
    return(NA_real_)
  }

  warn_labels <- c("partial", "best-effort", "best_effort", "legacy", "retired")
  fail_labels <- c("experimental", "unknown", "unsupported", "deprecated")

  if (any(labels %in% warn_labels)) {
    risk <- max(risk, 1)
  }
  if (any(labels %in% fail_labels)) {
    risk <- max(risk, 2)
  }

  risk
}

.library_card_read_count <- function(...) {
  frames <- list(...)
  read_ids <- unlist(lapply(frames, function(frame) {
    if (is.null(frame) || !is.data.frame(frame) || !"read_id" %in% names(frame)) {
      return(character())
    }
    as.character(frame$read_id)
  }), use.names = FALSE)
  read_ids <- unique(stats::na.omit(read_ids))
  if (length(read_ids) == 0) {
    return(NA_real_)
  }

  length(read_ids)
}

.library_card_row_matches <- function(data, pattern) {
  if (is.null(data) || nrow(data) == 0) {
    return(logical(0))
  }

  columns <- intersect(
    c("motif_kind", "motif_family", "motif_name", "motif_source"),
    names(data))
  if (length(columns) == 0) {
    return(rep(FALSE, nrow(data)))
  }

  Reduce(`|`, lapply(columns, function(column) {
    grepl(pattern, as.character(data[[column]]), ignore.case = TRUE)
  }))
}

.library_card_motif_names <- function(data) {
  if (is.null(data) || nrow(data) == 0) {
    return(character())
  }
  for (column in c("motif_name", "motif_family", "motif_kind")) {
    if (column %in% names(data)) {
      return(as.character(data[[column]]))
    }
  }

  character()
}

.library_card_as_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  label <- tolower(as.character(x))
  ifelse(label %in% c("true", "t", "yes", "1"), TRUE,
    ifelse(label %in% c("false", "f", "no", "0"), FALSE, NA))
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
