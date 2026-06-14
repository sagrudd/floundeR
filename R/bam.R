#' Summarise a BAM file for nanopore QC
#'
#' `bam_summary()` calls the curated `../bamana` Rust library in-process and
#' returns R-native tables for QC, review, provenance, and report assembly. The
#' binding is intentionally focused on read-only operational summary evidence:
#' it does not attempt to expose all Bamana commands.
#'
#' @param path Character scalar. Local BAM file to summarise.
#' @param sample_records Integer scalar. Use `0` for a full-file scan, or a
#'   positive value for a bounded evidence scan.
#' @param prefer_index Logical scalar. Prefer index-derived counts where Bamana
#'   can use a supported sidecar index.
#' @param include_mapq_hist Logical scalar. Include a MAPQ histogram table.
#' @param include_flags Logical scalar. Include flag-category evidence.
#' @param allow_incomplete Logical scalar. Allow Bamana to return partial
#'   evidence for incomplete BGZF/BAM input when supported by the backend.
#'
#' @return A named list with data frames including `status`, `evidence`,
#'   `header`, `counts`, `fractions`, `fractions_observed`, `mapq`, `mapping`,
#'   `anomalies`, `flag_categories`, `references`, `index_derived`, and
#'   `mapq_histogram`.
#'
#' @examples
#' \dontrun{
#' bam_summary("reads.bam")
#' }
#'
#' @export
bam_summary <- function(
  path,
  sample_records = 0L,
  prefer_index = FALSE,
  include_mapq_hist = FALSE,
  include_flags = TRUE,
  allow_incomplete = FALSE
) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_bam_error("`path` must be a non-missing character scalar.",
      category = "path"
    )
  }
  sample_records <- .flounder_nonnegative_integer_scalar(
    sample_records,
    "sample_records"
  )
  prefer_index <- .flounder_logical_scalar(prefer_index, "prefer_index")
  include_mapq_hist <- .flounder_logical_scalar(
    include_mapq_hist,
    "include_mapq_hist"
  )
  include_flags <- .flounder_logical_scalar(include_flags, "include_flags")
  allow_incomplete <- .flounder_logical_scalar(
    allow_incomplete,
    "allow_incomplete"
  )

  response <- .flounder_bam_summary_call(
    path,
    sample_records,
    prefer_index,
    include_mapq_hist,
    include_flags,
    allow_incomplete
  )
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_bam_error(
      response$error %||% "BAM summary failed in the Rust extension.",
      category = response$category %||% "unknown",
      code = response$code %||% NA_character_,
      detail = .flounder_empty_scalar_to_na(response$detail %||% NA_character_),
      hint = .flounder_empty_scalar_to_na(response$hint %||% NA_character_)
    )
  }

  .flounder_bam_summary_normalise(response$data)
}

#' Verify a BAM file header and container identity
#'
#' `bam_verify()` calls Bamana's shallow verification library path in-process.
#' It confirms BGZF container recognition, BAM magic, and native BAM header
#' parsing, but it does not scan alignment records and does not prove EOF
#' completeness.
#'
#' @param path Character scalar. Local BAM file to verify.
#'
#' @return A one-row data frame with command status, detected format,
#'   container, shallow/deep validation flags, checks performed, confidence,
#'   semantic note, and Bamana error metadata columns.
#'
#' @examples
#' \dontrun{
#' bam_verify("reads.bam")
#' }
#'
#' @export
bam_verify <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_bam_error("`path` must be a non-missing character scalar.",
      category = "path"
    )
  }

  response <- .flounder_bam_verify_call(path)
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_bam_error(
      response$error %||% "BAM verification failed in the Rust extension.",
      category = response$category %||% "unknown",
      code = response$code %||% NA_character_,
      detail = .flounder_empty_scalar_to_na(response$detail %||% NA_character_),
      hint = .flounder_empty_scalar_to_na(response$hint %||% NA_character_)
    )
  }

  .flounder_bam_command_table_normalise(response$data,
    class = "floundeR_bam_verify"
  )
}

#' Validate BAM structure and selected consistency checks
#'
#' `bam_validate()` calls Bamana's structural validation library path
#' in-process. Validation findings are returned as evidence tables even when
#' Bamana reports validation failure; missing or unreadable inputs without a
#' validation payload still raise typed Bamana-backed R conditions.
#'
#' @param path Character scalar. Local BAM file to validate.
#' @param max_errors Positive integer scalar. Maximum error findings to retain.
#' @param max_warnings Non-negative integer scalar. Maximum warning findings to
#'   retain.
#' @param header_only Logical scalar. Validate the header only.
#' @param records Optional positive integer scalar. Validate at most this many
#'   records after the header.
#' @param fail_fast Logical scalar. Stop at the first error-level finding.
#' @param include_warnings Logical scalar. Include warning-level findings.
#'
#' @return A named list with `status`, `summary`, `findings`, and `error` data
#'   frames. `status$ok` preserves Bamana's command envelope result; a
#'   validation failure with findings therefore returns `ok = FALSE` rather
#'   than raising an R error.
#'
#' @examples
#' \dontrun{
#' bam_validate("reads.bam")
#' }
#'
#' @export
bam_validate <- function(
  path,
  max_errors = 100L,
  max_warnings = 100L,
  header_only = FALSE,
  records = NULL,
  fail_fast = FALSE,
  include_warnings = TRUE
) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_bam_error("`path` must be a non-missing character scalar.",
      category = "path"
    )
  }
  max_errors <- .flounder_positive_bam_integer_scalar(max_errors, "max_errors")
  max_warnings <- .flounder_nonnegative_integer_scalar(
    max_warnings,
    "max_warnings"
  )
  header_only <- .flounder_logical_scalar(header_only, "header_only")
  records <- .flounder_optional_positive_bam_integer_scalar(records, "records")
  fail_fast <- .flounder_logical_scalar(fail_fast, "fail_fast")
  include_warnings <- .flounder_logical_scalar(
    include_warnings,
    "include_warnings"
  )

  response <- .flounder_bam_validate_call(
    path,
    max_errors,
    max_warnings,
    header_only,
    records,
    fail_fast,
    include_warnings
  )
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_bam_error(
      response$error %||% "BAM validation failed in the Rust extension.",
      category = response$category %||% "unknown",
      code = response$code %||% NA_character_,
      detail = .flounder_empty_scalar_to_na(response$detail %||% NA_character_),
      hint = .flounder_empty_scalar_to_na(response$hint %||% NA_character_)
    )
  }

  .flounder_bam_validate_normalise(response$data)
}

#' Check canonical BGZF EOF evidence
#'
#' `bam_check_eof()` reports canonical BGZF EOF-marker evidence through
#' Bamana's in-process Rust library path. It reports tail completeness only and
#' does not imply BAM header, record, or aux-field validity.
#'
#' @param path Character scalar. Local BAM file to inspect.
#'
#' @return A one-row data frame with command status, detected format,
#'   `bgzf_eof_present`, `complete`, semantic note, and Bamana error metadata
#'   columns. Missing EOF is returned as `complete = FALSE` evidence.
#'
#' @examples
#' \dontrun{
#' bam_check_eof("reads.bam")
#' }
#'
#' @export
bam_check_eof <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_bam_error("`path` must be a non-missing character scalar.",
      category = "path"
    )
  }

  response <- .flounder_bam_check_eof_call(path)
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_bam_error(
      response$error %||% "BAM EOF check failed in the Rust extension.",
      category = response$category %||% "unknown",
      code = response$code %||% NA_character_,
      detail = .flounder_empty_scalar_to_na(response$detail %||% NA_character_),
      hint = .flounder_empty_scalar_to_na(response$hint %||% NA_character_)
    )
  }

  .flounder_bam_command_table_normalise(response$data,
    class = "floundeR_bam_check_eof"
  )
}

#' Check BAM index evidence for QC report cards
#'
#' `bam_check_index()` calls Bamana's in-process Rust index inspection path and
#' returns stable R tables for index presence, support level, freshness, and
#' candidate sidecars. It is read-only and does not create or repair indexes.
#'
#' @param path Character scalar. Local BAM file to inspect.
#' @param require Logical scalar. Treat a missing index as a Bamana command
#'   failure while still returning structured index evidence when available.
#' @param prefer_csi Logical scalar. Prefer CSI over BAI when both sidecars are
#'   present and usable.
#'
#' @return A named list with `status`, `index`, `candidates`, and `error` data
#'   frames.
#'
#' @examples
#' \dontrun{
#' bam_check_index("reads.bam")
#' }
#'
#' @export
bam_check_index <- function(path, require = FALSE, prefer_csi = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_bam_error("`path` must be a non-missing character scalar.",
      category = "path"
    )
  }
  require <- .flounder_logical_scalar(require, "require")
  prefer_csi <- .flounder_logical_scalar(prefer_csi, "prefer_csi")

  response <- .flounder_bam_check_index_call(path, require, prefer_csi)
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_bam_error(
      response$error %||% "BAM index check failed in the Rust extension.",
      category = response$category %||% "unknown",
      code = response$code %||% NA_character_,
      detail = .flounder_empty_scalar_to_na(response$detail %||% NA_character_),
      hint = .flounder_empty_scalar_to_na(response$hint %||% NA_character_)
    )
  }

  .flounder_bam_check_index_normalise(response$data)
}

#' Check whether BAM mapping evidence is present
#'
#' `bam_check_map()` calls Bamana's in-process Rust mapping check. It reports
#' whether mapped reads are evidenced by a supported index or by a bounded/full
#' record scan, with reference-level evidence for QC review.
#'
#' @param path Character scalar. Local BAM file to inspect.
#' @param sample_records Integer scalar. Use `0` for a full-file scan, or a
#'   positive value for a bounded evidence scan.
#' @param prefer_index Logical scalar. Prefer supported index-derived evidence
#'   where possible.
#'
#' @return A named list with `status`, `index`, `summary`, and `references`
#'   data frames.
#'
#' @examples
#' \dontrun{
#' bam_check_map("reads.bam", sample_records = 1000)
#' }
#'
#' @export
bam_check_map <- function(path, sample_records = 1000L, prefer_index = TRUE) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_bam_error("`path` must be a non-missing character scalar.",
      category = "path"
    )
  }
  sample_records <- .flounder_nonnegative_integer_scalar(
    sample_records,
    "sample_records"
  )
  prefer_index <- .flounder_logical_scalar(prefer_index, "prefer_index")

  response <- .flounder_bam_check_map_call(path, sample_records, prefer_index)
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_bam_error(
      response$error %||% "BAM mapping check failed in the Rust extension.",
      category = response$category %||% "unknown",
      code = response$code %||% NA_character_,
      detail = .flounder_empty_scalar_to_na(response$detail %||% NA_character_),
      hint = .flounder_empty_scalar_to_na(response$hint %||% NA_character_)
    )
  }

  .flounder_bam_check_map_normalise(response$data)
}

#' Check BAM sort declaration and observed sort evidence
#'
#' `bam_check_sort()` calls Bamana's in-process Rust sort check and returns a
#' one-row QC table comparing declared header order with observed record order.
#'
#' @param path Character scalar. Local BAM file to inspect.
#' @param sample_records Integer scalar. Use `0` for a full-file scan, or a
#'   positive value for a bounded evidence scan.
#' @param strict Logical scalar. Ask Bamana to use strict sort interpretation.
#'
#' @return A one-row data frame with declared sort order, observed evidence,
#'   agreement, confidence, and first-violation columns.
#'
#' @examples
#' \dontrun{
#' bam_check_sort("reads.bam")
#' }
#'
#' @export
bam_check_sort <- function(path, sample_records = 1000L, strict = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_bam_error("`path` must be a non-missing character scalar.",
      category = "path"
    )
  }
  sample_records <- .flounder_nonnegative_integer_scalar(
    sample_records,
    "sample_records"
  )
  strict <- .flounder_logical_scalar(strict, "strict")

  response <- .flounder_bam_check_sort_call(path, sample_records, strict)
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_bam_error(
      response$error %||% "BAM sort check failed in the Rust extension.",
      category = response$category %||% "unknown",
      code = response$code %||% NA_character_,
      detail = .flounder_empty_scalar_to_na(response$detail %||% NA_character_),
      hint = .flounder_empty_scalar_to_na(response$hint %||% NA_character_)
    )
  }

  .flounder_bam_check_sort_normalise(response$data)
}

#' Check BAM aux-tag evidence for library/QC review
#'
#' `bam_check_tag()` calls Bamana's in-process Rust aux-tag check. It is useful
#' for QC report cards that need bounded or full-scan evidence for tags such as
#' read groups, molecular identifiers, or alignment annotations.
#'
#' @param path Character scalar. Local BAM file to inspect.
#' @param tag Character scalar. Two-character BAM aux tag to search for.
#' @param sample_records Integer scalar. Number of records to examine when
#'   `full_scan = FALSE`; use `0` only with `full_scan = TRUE`.
#' @param full_scan Logical scalar. Scan the full file for absence/presence
#'   evidence.
#' @param require_type Optional character scalar. Bamana aux type code to
#'   require, for example `"Z"` or `"i"`.
#' @param count_hits Logical scalar. Count all matching records during the scan.
#'
#' @return A one-row data frame with tag evidence, mode, confidence, and Bamana
#'   error metadata columns.
#'
#' @examples
#' \dontrun{
#' bam_check_tag("reads.bam", "RG", full_scan = TRUE)
#' }
#'
#' @export
bam_check_tag <- function(
  path,
  tag,
  sample_records = 1000L,
  full_scan = FALSE,
  require_type = NULL,
  count_hits = FALSE
) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_bam_error("`path` must be a non-missing character scalar.",
      category = "path"
    )
  }
  if (!is.character(tag) || length(tag) != 1L || is.na(tag)) {
    .flounder_bam_error("`tag` must be a non-missing character scalar.",
      category = "argument"
    )
  }
  sample_records <- .flounder_nonnegative_integer_scalar(
    sample_records,
    "sample_records"
  )
  full_scan <- .flounder_logical_scalar(full_scan, "full_scan")
  require_type <- .flounder_optional_character_scalar(require_type, "require_type")
  count_hits <- .flounder_logical_scalar(count_hits, "count_hits")

  response <- .flounder_bam_check_tag_call(
    path,
    tag,
    sample_records,
    full_scan,
    require_type,
    count_hits
  )
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_bam_error(
      response$error %||% "BAM tag check failed in the Rust extension.",
      category = response$category %||% "unknown",
      code = response$code %||% NA_character_,
      detail = .flounder_empty_scalar_to_na(response$detail %||% NA_character_),
      hint = .flounder_empty_scalar_to_na(response$hint %||% NA_character_)
    )
  }

  .flounder_bam_check_tag_normalise(response$data)
}

.flounder_bam_summary_call <- function(
  path,
  sample_records,
  prefer_index,
  include_mapq_hist,
  include_flags,
  allow_incomplete
) {
  .Call(
    "flounder_bam_summary",
    path,
    as.integer(sample_records),
    prefer_index,
    include_mapq_hist,
    include_flags,
    allow_incomplete,
    PACKAGE = "floundeR"
  )
}

.flounder_bam_verify_call <- function(path) {
  .Call("flounder_bam_verify", path, PACKAGE = "floundeR")
}

.flounder_bam_validate_call <- function(
  path,
  max_errors,
  max_warnings,
  header_only,
  records,
  fail_fast,
  include_warnings
) {
  .Call(
    "flounder_bam_validate",
    path,
    as.integer(max_errors),
    as.integer(max_warnings),
    header_only,
    .flounder_optional_integer_label(records),
    fail_fast,
    include_warnings,
    PACKAGE = "floundeR"
  )
}

.flounder_bam_check_eof_call <- function(path) {
  .Call("flounder_bam_check_eof", path, PACKAGE = "floundeR")
}

.flounder_bam_check_index_call <- function(path, require, prefer_csi) {
  .Call(
    "flounder_bam_check_index",
    path,
    require,
    prefer_csi,
    PACKAGE = "floundeR"
  )
}

.flounder_bam_check_map_call <- function(path, sample_records, prefer_index) {
  .Call(
    "flounder_bam_check_map",
    path,
    as.integer(sample_records),
    prefer_index,
    PACKAGE = "floundeR"
  )
}

.flounder_bam_check_sort_call <- function(path, sample_records, strict) {
  .Call(
    "flounder_bam_check_sort",
    path,
    as.integer(sample_records),
    strict,
    PACKAGE = "floundeR"
  )
}

.flounder_bam_check_tag_call <- function(
  path,
  tag,
  sample_records,
  full_scan,
  require_type,
  count_hits
) {
  .Call(
    "flounder_bam_check_tag",
    path,
    tag,
    as.integer(sample_records),
    full_scan,
    require_type %||% "",
    count_hits,
    PACKAGE = "floundeR"
  )
}

.flounder_bam_summary_normalise <- function(data) {
  if (!is.list(data)) {
    .flounder_bam_error("BAM summary returned an unsupported payload shape.")
  }

  for (name in names(data)) {
    data[[name]] <- as.data.frame(data[[name]], stringsAsFactors = FALSE)
  }

  data$status$schema_version <- as.integer(data$status$schema_version)
  data$status$analysis_wall_seconds <- .flounder_nan_to_na(
    data$status$analysis_wall_seconds
  )
  data$status$confidence <- .flounder_empty_to_na(data$status$confidence)
  data$status$semantic_note <- .flounder_empty_to_na(data$status$semantic_note)

  data$evidence$schema_version <- as.integer(data$evidence$schema_version)
  data$evidence$records_scanned <- as.numeric(data$evidence$records_scanned)

  data$header$schema_version <- as.integer(data$header$schema_version)
  data$header$references_defined <- as.integer(data$header$references_defined)
  for (column in c("sort_order", "sub_sort_order", "group_order")) {
    data$header[[column]] <- .flounder_empty_to_na(data$header[[column]])
  }

  data$counts$schema_version <- as.integer(data$counts$schema_version)
  for (column in setdiff(names(data$counts), "schema_version")) {
    data$counts[[column]] <- .flounder_nan_to_na(data$counts[[column]])
  }

  data$fractions <- .flounder_bam_fraction_table(data$fractions)
  data$fractions_observed <- .flounder_bam_fraction_table(
    data$fractions_observed
  )

  data$mapq$schema_version <- as.integer(data$mapq$schema_version)
  for (column in setdiff(names(data$mapq), "schema_version")) {
    data$mapq[[column]] <- .flounder_nan_to_na(data$mapq[[column]])
  }

  data$mapping$schema_version <- as.integer(data$mapping$schema_version)
  for (column in c(
    "references_with_mapped_reads",
    "references_with_mapped_reads_observed"
  )) {
    data$mapping[[column]] <- .flounder_nan_to_na(data$mapping[[column]])
  }

  data$anomalies$schema_version <- as.integer(data$anomalies$schema_version)
  data$anomalies$contradictory_mapping_state_records <- .flounder_nan_to_na(
    data$anomalies$contradictory_mapping_state_records
  )

  data$flag_categories$schema_version <- as.integer(
    data$flag_categories$schema_version
  )
  for (column in setdiff(names(data$flag_categories), "schema_version")) {
    data$flag_categories[[column]] <- .flounder_nan_to_na(
      data$flag_categories[[column]]
    )
  }

  data$references$schema_version <- as.integer(data$references$schema_version)
  for (column in c("length", "mapped_reads", "unmapped_reads")) {
    data$references[[column]] <- .flounder_nan_to_na(data$references[[column]])
  }
  data$references$observed_mapped <- .flounder_empty_to_na(
    data$references$observed_mapped
  )

  data$index_derived$schema_version <- as.integer(
    data$index_derived$schema_version
  )
  for (column in c(
    "kind",
    "total_mapped_reads",
    "total_unmapped_reads",
    "references_with_mapped_reads",
    "note"
  )) {
    if (is.character(data$index_derived[[column]])) {
      data$index_derived[[column]] <- .flounder_empty_to_na(
        data$index_derived[[column]]
      )
    } else {
      data$index_derived[[column]] <- .flounder_nan_to_na(
        data$index_derived[[column]]
      )
    }
  }

  data$mapq_histogram$schema_version <- as.integer(
    data$mapq_histogram$schema_version
  )
  data$mapq_histogram$mapq <- as.integer(data$mapq_histogram$mapq)
  data$mapq_histogram$read_count <- .flounder_nan_to_na(
    data$mapq_histogram$read_count
  )

  class(data) <- c("floundeR_bam_summary", "list")
  data
}

.flounder_bam_command_table_normalise <- function(data, class) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  data$schema_version <- as.integer(data$schema_version)
  data$analysis_wall_seconds <- .flounder_nan_to_na(
    data$analysis_wall_seconds
  )
  for (column in c(
    "checks_performed",
    "semantic_note",
    "error_code",
    "error_message",
    "error_detail",
    "error_hint"
  )) {
    if (column %in% names(data)) {
      data[[column]] <- .flounder_empty_to_na(data[[column]])
    }
  }
  class(data) <- c(class, "data.frame")
  data
}

.flounder_bam_validate_normalise <- function(data) {
  if (!is.list(data)) {
    .flounder_bam_error("BAM validation returned an unsupported payload shape.")
  }

  for (name in names(data)) {
    data[[name]] <- as.data.frame(data[[name]], stringsAsFactors = FALSE)
  }

  data$status$schema_version <- as.integer(data$status$schema_version)
  data$status$analysis_wall_seconds <- .flounder_nan_to_na(
    data$status$analysis_wall_seconds
  )
  data$status$semantic_note <- .flounder_empty_to_na(data$status$semantic_note)

  data$summary$schema_version <- as.integer(data$summary$schema_version)
  for (column in c("records_examined", "errors", "warnings", "infos")) {
    data$summary[[column]] <- .flounder_nan_to_na(data$summary[[column]])
  }

  data$findings$schema_version <- as.integer(data$findings$schema_version)
  data$findings$record_index <- .flounder_nan_to_na(data$findings$record_index)
  data$findings$reference_name <- .flounder_empty_to_na(
    data$findings$reference_name
  )
  data$findings$tag <- .flounder_empty_to_na(data$findings$tag)

  data$error$schema_version <- as.integer(data$error$schema_version)
  for (column in c("code", "message", "detail", "hint")) {
    data$error[[column]] <- .flounder_empty_to_na(data$error[[column]])
  }

  class(data) <- c("floundeR_bam_validate", "list")
  data
}

.flounder_bam_check_index_normalise <- function(data) {
  if (!is.list(data)) {
    .flounder_bam_error("BAM index check returned an unsupported payload shape.")
  }

  for (name in names(data)) {
    data[[name]] <- as.data.frame(data[[name]], stringsAsFactors = FALSE)
  }

  data$status <- .flounder_bam_status_table_normalise(data$status)

  data$index$schema_version <- as.integer(data$index$schema_version)
  for (column in c(
    "selected_path",
    "kind",
    "syntactically_valid",
    "stale",
    "bam_newer_than_index",
    "notes"
  )) {
    data$index[[column]] <- .flounder_empty_to_na(data$index[[column]])
  }

  data$candidates$schema_version <- as.integer(data$candidates$schema_version)

  data$error$schema_version <- as.integer(data$error$schema_version)
  for (column in c("code", "message", "detail", "hint")) {
    data$error[[column]] <- .flounder_empty_to_na(data$error[[column]])
  }

  class(data) <- c("floundeR_bam_check_index", "list")
  data
}

.flounder_bam_check_map_normalise <- function(data) {
  if (!is.list(data)) {
    .flounder_bam_error("BAM mapping check returned an unsupported payload shape.")
  }

  for (name in names(data)) {
    data[[name]] <- as.data.frame(data[[name]], stringsAsFactors = FALSE)
  }

  data$status <- .flounder_bam_status_table_normalise(data$status)
  data$status$has_mapped_reads <- .flounder_empty_to_na(
    data$status$has_mapped_reads
  )

  data$index$schema_version <- as.integer(data$index$schema_version)
  data$index$kind <- .flounder_empty_to_na(data$index$kind)
  data$index$diagnostic_detail <- .flounder_empty_to_na(
    data$index$diagnostic_detail
  )

  data$summary$schema_version <- as.integer(data$summary$schema_version)
  for (column in setdiff(names(data$summary), "schema_version")) {
    data$summary[[column]] <- .flounder_nan_to_na(data$summary[[column]])
  }

  data$references$schema_version <- as.integer(data$references$schema_version)
  for (column in c("length", "mapped_reads", "unmapped_reads")) {
    data$references[[column]] <- .flounder_nan_to_na(data$references[[column]])
  }
  data$references$observed <- .flounder_empty_to_na(data$references$observed)

  class(data) <- c("floundeR_bam_check_map", "list")
  data
}

.flounder_bam_check_sort_normalise <- function(data) {
  data <- .flounder_bam_command_table_normalise(data,
    class = "floundeR_bam_check_sort"
  )
  for (column in c(
    "declared_so",
    "declared_ss",
    "declared_go",
    "observed_sub_order",
    "appears_sorted",
    "first_violation_reason",
    "header_matches_observation"
  )) {
    data[[column]] <- .flounder_empty_to_na(data[[column]])
  }
  data$records_examined <- .flounder_nan_to_na(data$records_examined)
  data$first_violation_record_index <- .flounder_nan_to_na(
    data$first_violation_record_index
  )
  data
}

.flounder_bam_check_tag_normalise <- function(data) {
  data <- .flounder_bam_command_table_normalise(data,
    class = "floundeR_bam_check_tag"
  )
  for (column in c("required_type", "confidence", "semantic_note")) {
    data[[column]] <- .flounder_empty_to_na(data[[column]])
  }
  data$records_examined <- .flounder_nan_to_na(data$records_examined)
  data$records_with_tag <- .flounder_nan_to_na(data$records_with_tag)
  data
}

.flounder_bam_status_table_normalise <- function(data) {
  data$schema_version <- as.integer(data$schema_version)
  data$analysis_wall_seconds <- .flounder_nan_to_na(
    data$analysis_wall_seconds
  )
  for (column in c("semantic_note", "confidence")) {
    if (column %in% names(data)) {
      data[[column]] <- .flounder_empty_to_na(data[[column]])
    }
  }
  data
}

.flounder_bam_fraction_table <- function(data) {
  data$schema_version <- as.integer(data$schema_version)
  for (column in setdiff(names(data), c("schema_version", "scope"))) {
    data[[column]] <- .flounder_nan_to_na(data[[column]])
  }
  data
}

.flounder_nonnegative_integer_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0L || x != as.integer(x)) {
    .flounder_bam_error(
      paste0("`", name, "` must be zero or a positive integer scalar."),
      category = "argument"
    )
  }
  as.integer(x)
}

.flounder_positive_bam_integer_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1L || x != as.integer(x)) {
    .flounder_bam_error(
      paste0("`", name, "` must be a positive integer scalar."),
      category = "argument"
    )
  }
  as.integer(x)
}

.flounder_optional_positive_bam_integer_scalar <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }
  .flounder_positive_bam_integer_scalar(x, name)
}

.flounder_optional_character_scalar <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    .flounder_bam_error(
      paste0("`", name, "` must be NULL or a non-missing character scalar."),
      category = "argument"
    )
  }
  x
}

.flounder_logical_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .flounder_bam_error(
      paste0("`", name, "` must be TRUE or FALSE."),
      category = "argument"
    )
  }
  x
}

.flounder_empty_scalar_to_na <- function(x) {
  if (length(x) != 1L || is.na(x) || identical(x, "")) {
    return(NA_character_)
  }
  x
}

.flounder_bam_error <- function(
  message,
  category = "unknown",
  code = NA_character_,
  detail = NA_character_,
  hint = NA_character_,
  parent = NULL,
  call = NULL
) {
  subclass <- switch(
    category,
    path = "floundeR_bam_path_error",
    format = "floundeR_bam_format_error",
    index = "floundeR_bam_index_error",
    argument = "floundeR_bam_argument_error",
    "floundeR_bam_unknown_error"
  )
  condition <- structure(
    list(
      message = message,
      code = code,
      detail = detail,
      hint = hint,
      parent = parent,
      call = call
    ),
    class = c(
      subclass,
      "floundeR_bam_error",
      "floundeR_error",
      "error",
      "condition"
    )
  )
  stop(condition)
}
