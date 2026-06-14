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
