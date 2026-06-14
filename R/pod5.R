#' Discover local POD5-containing folders
#'
#' `pod5_find()` recursively searches a local directory tree for folders that
#' directly contain one or more `.pod5` files. Discovery is implemented through
#' the in-process Rust extension and the curated `pod5-tools` library API; it
#' does not invoke an external command-line process.
#'
#' @param path Character scalar. Directory to search.
#'
#' @return A data frame with one row per discovered directory and the columns
#'   `path`, `pod5_file_count`, `total_bytes`, `oldest_modified_utc`, and
#'   `newest_modified_utc`.
#'
#' @examples
#' run_dir <- file.path(tempdir(), "flounder-pod5-find-example")
#' unlink(run_dir, recursive = TRUE)
#' dir.create(file.path(run_dir, "sample-a"), recursive = TRUE)
#' dir.create(file.path(run_dir, "sample-b"), recursive = TRUE)
#'
#' writeBin(charToRaw("pod5 placeholder"), file.path(run_dir, "sample-a", "a.pod5"))
#' writeBin(charToRaw("ignored"), file.path(run_dir, "sample-a", "a.fastq"))
#' writeBin(charToRaw("pod5 placeholder"), file.path(run_dir, "sample-b", "b.pod5"))
#'
#' pod5_find(run_dir)
#'
#' @export
pod5_find <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_pod5_error("`path` must be a non-missing character scalar.")
  }

  response <- .flounder_pod5_find_call(path)
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_pod5_error(
      response$error %||% "POD5 discovery failed in the Rust extension."
    )
  }

  data <- response$data
  names(data) <- c(
    "path",
    "pod5_file_count",
    "total_bytes",
    "oldest_modified_utc",
    "newest_modified_utc"
  )
  data$pod5_file_count <- as.integer(data$pod5_file_count)
  data$total_bytes <- as.numeric(data$total_bytes)
  data
}

#' Verify a local POD5 candidate file
#'
#' `pod5_verify()` checks a single local file through the in-process Rust
#' extension and the curated `pod5-tools` library API. The current backend
#' performs fast extension and fixed POD5 signature checks, then reports deeper
#' layout/schema checks as `not_checked` until the full parser backend is
#' connected. A signature-valid file therefore returns `overall_status =
#' "incomplete"` rather than overstating full POD5 validity.
#'
#' @param path Character scalar. Candidate POD5 file to verify.
#'
#' @return A data frame with one row per verification check and the columns
#'   `path`, `size_bytes`, `overall_status`, `check`, `category`, `status`, and
#'   `detail`.
#'
#' @examples
#' signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
#' path <- file.path(tempdir(), "flounder-signature-example.pod5")
#' writeBin(c(signature, as.raw(rep(0, 16)), signature), path)
#' pod5_verify(path)
#'
#' @export
pod5_verify <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_pod5_error("`path` must be a non-missing character scalar.",
      category = "path"
    )
  }

  response <- .flounder_pod5_verify_call(path)
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_pod5_error(
      response$error %||% "POD5 verification failed in the Rust extension.",
      category = response$category %||% "unknown"
    )
  }

  data <- response$data
  names(data) <- c(
    "path",
    "size_bytes",
    "overall_status",
    "check",
    "category",
    "status",
    "detail"
  )
  data$size_bytes <- as.numeric(data$size_bytes)
  data
}

#' Inspect local POD5 file metadata
#'
#' `pod5_file_info()` returns file-level POD5 metadata through the in-process
#' Rust extension and the curated `pod5-tools` library API. Until the deep POD5
#' parser backend is connected, internal fields such as read count, flow cell
#' ID, sequencing kit, and POD5 version are returned as missing values, while
#' file size and integrity availability are still reported.
#'
#' @param path Character scalar. POD5 file to inspect.
#'
#' @return A one-row data frame with the columns `path`, `size_bytes`,
#'   `flow_cell_id`, `sequencing_kit`, `read_count`, `acquisition_start_utc`,
#'   `duration_seconds`, `pod5_version`, `integrity_status`, and
#'   `integrity_reason`.
#'
#' @examples
#' signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
#' path <- file.path(tempdir(), "flounder-file-info-example.pod5")
#' writeBin(c(signature, as.raw(rep(0, 16)), signature), path)
#' pod5_file_info(path)
#'
#' @export
pod5_file_info <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_pod5_error("`path` must be a non-missing character scalar.",
      category = "path"
    )
  }

  response <- .flounder_pod5_file_info_call(path)
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_pod5_error(
      response$error %||% "POD5 file metadata inspection failed in the Rust extension.",
      category = response$category %||% "unknown"
    )
  }

  data <- response$data
  names(data) <- c(
    "path",
    "size_bytes",
    "flow_cell_id",
    "sequencing_kit",
    "read_count",
    "acquisition_start_utc",
    "duration_seconds",
    "pod5_version",
    "integrity_status",
    "integrity_reason"
  )
  data$size_bytes <- as.numeric(data$size_bytes)
  data$read_count <- .flounder_nan_to_na(data$read_count)
  data$duration_seconds <- .flounder_nan_to_na(data$duration_seconds)
  for (column in c(
    "flow_cell_id",
    "sequencing_kit",
    "acquisition_start_utc",
    "pod5_version",
    "integrity_reason"
  )) {
    data[[column]] <- .flounder_empty_to_na(data[[column]])
  }
  data
}

#' Summarise a local POD5 folder or run tree
#'
#' `pod5_folder_info()` aggregates file-level POD5 metadata and fast
#' verification state across a local directory tree through the in-process Rust
#' extension and the curated `pod5-tools` library API. The current backend
#' detects duplicate file names and implemented verification failures while
#' reporting unavailable deep POD5 metadata explicitly in the `warnings` and
#' `integrity_*` fields.
#'
#' @param path Character scalar. Directory containing POD5 files.
#'
#' @return A one-row data frame with the columns `path`, `pod5_file_count`,
#'   `total_bytes`, `total_reads`, `flow_cell_ids`, `sequencing_kits`,
#'   `acquisition_start_utc`, `acquisition_end_utc`, `integrity_status`,
#'   `integrity_reason`, `failed_file_count`, `verification_failed_count`,
#'   `duplicate_file_names`, and `warnings`.
#'
#' @examples
#' signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
#' run_dir <- file.path(tempdir(), "flounder-folder-info-example")
#' unlink(run_dir, recursive = TRUE)
#' dir.create(run_dir, recursive = TRUE)
#' writeBin(c(signature, as.raw(rep(0, 16)), signature), file.path(run_dir, "a.pod5"))
#' pod5_folder_info(run_dir)
#'
#' @export
pod5_folder_info <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_pod5_error("`path` must be a non-missing character scalar.",
      category = "path"
    )
  }

  response <- .flounder_pod5_folder_info_call(path)
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_pod5_error(
      response$error %||% "POD5 folder inspection failed in the Rust extension.",
      category = response$category %||% "unknown"
    )
  }

  data <- response$data
  names(data) <- c(
    "path",
    "pod5_file_count",
    "total_bytes",
    "total_reads",
    "flow_cell_ids",
    "sequencing_kits",
    "acquisition_start_utc",
    "acquisition_end_utc",
    "integrity_status",
    "integrity_reason",
    "failed_file_count",
    "verification_failed_count",
    "duplicate_file_names",
    "warnings"
  )
  data$pod5_file_count <- as.integer(data$pod5_file_count)
  data$total_bytes <- as.numeric(data$total_bytes)
  data$total_reads <- .flounder_nan_to_na(data$total_reads)
  data$failed_file_count <- as.integer(data$failed_file_count)
  data$verification_failed_count <- as.integer(data$verification_failed_count)
  for (column in c(
    "flow_cell_ids",
    "sequencing_kits",
    "acquisition_start_utc",
    "acquisition_end_utc",
    "integrity_reason",
    "duplicate_file_names",
    "warnings"
  )) {
    data[[column]] <- .flounder_empty_to_na(data[[column]])
  }
  data
}

#' Build a versioned POD5 collection manifest
#'
#' `pod5_manifest()` inventories a local POD5 file, folder, or run tree through
#' the in-process Rust extension and the curated `pod5-tools` library API. The
#' returned manifest is intentionally row-oriented for R and synoptikon handoff:
#' every row carries the manifest schema version and source path alongside one
#' POD5 entry.
#'
#' @param path Character scalar. POD5 file or directory to inventory.
#'
#' @return A data frame with one row per POD5 file and the columns
#'   `schema_version`, `source`, `relative_path`, `path`, `size_bytes`,
#'   `verification_status`, and `verification_failed_checks`.
#'
#' @examples
#' signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
#' run_dir <- file.path(tempdir(), "flounder-manifest-example")
#' unlink(run_dir, recursive = TRUE)
#' dir.create(run_dir, recursive = TRUE)
#' writeBin(c(signature, as.raw(rep(0, 16)), signature), file.path(run_dir, "a.pod5"))
#' pod5_manifest(run_dir)
#'
#' @export
pod5_manifest <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_pod5_error("`path` must be a non-missing character scalar.",
      category = "path"
    )
  }

  response <- .flounder_pod5_manifest_call(path)
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_pod5_error(
      response$error %||% "POD5 manifest construction failed in the Rust extension.",
      category = response$category %||% "unknown"
    )
  }

  data <- response$data
  names(data) <- c(
    "schema_version",
    "source",
    "relative_path",
    "path",
    "size_bytes",
    "verification_status",
    "verification_failed_checks"
  )
  data$schema_version <- as.integer(data$schema_version)
  data$size_bytes <- as.numeric(data$size_bytes)
  data$verification_failed_checks <- as.integer(data$verification_failed_checks)
  data
}

#' Compare two POD5 collections or manifests
#'
#' `pod5_compare()` compares two local POD5 files, folders, run trees, or
#' `pod5-tools` manifest JSON files through the in-process Rust extension. It is
#' designed for operational handoff checks: missing files and changed file-size
#' or verification-state differences are returned as rows rather than rendered
#' text.
#'
#' @param left,right Character scalars. POD5 inputs or manifest JSON files to
#'   compare.
#'
#' @return A data frame with the columns `status`, `kind`, `relative_path`,
#'   `left_size_bytes`, `right_size_bytes`, `left_verification_status`, and
#'   `right_verification_status`. Matching inputs return one row with
#'   `status = "match"` and `kind = "match"`.
#'
#' @examples
#' signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
#' left <- file.path(tempdir(), "flounder-compare-left")
#' right <- file.path(tempdir(), "flounder-compare-right")
#' unlink(c(left, right), recursive = TRUE)
#' dir.create(left, recursive = TRUE)
#' dir.create(right, recursive = TRUE)
#' writeBin(c(signature, as.raw(rep(0, 16)), signature), file.path(left, "a.pod5"))
#' writeBin(c(signature, as.raw(rep(1, 16)), signature), file.path(right, "a.pod5"))
#' pod5_compare(left, right)
#'
#' @export
pod5_compare <- function(left, right) {
  if (!is.character(left) || length(left) != 1L || is.na(left)) {
    .flounder_pod5_error("`left` must be a non-missing character scalar.",
      category = "path"
    )
  }
  if (!is.character(right) || length(right) != 1L || is.na(right)) {
    .flounder_pod5_error("`right` must be a non-missing character scalar.",
      category = "path"
    )
  }

  response <- .flounder_pod5_compare_call(left, right)
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_pod5_error(
      response$error %||% "POD5 comparison failed in the Rust extension.",
      category = response$category %||% "unknown"
    )
  }

  data <- response$data
  names(data) <- c(
    "status",
    "kind",
    "relative_path",
    "left_size_bytes",
    "right_size_bytes",
    "left_verification_status",
    "right_verification_status"
  )
  data$relative_path <- .flounder_empty_to_na(data$relative_path)
  data$left_size_bytes <- .flounder_nan_to_na(data$left_size_bytes)
  data$right_size_bytes <- .flounder_nan_to_na(data$right_size_bytes)
  data$left_verification_status <- .flounder_empty_to_na(data$left_verification_status)
  data$right_verification_status <- .flounder_empty_to_na(data$right_verification_status)
  data
}

#' Plan a read-only POD5 collection subdivision
#'
#' `pod5_subdivide_plan()` asks the in-process Rust extension and curated
#' `pod5-tools` library API to build a deterministic, read-only plan describing
#' how a POD5 file, folder, run tree, or manifest could be split for
#' demonstration, review, or workflow development. This function does not write
#' POD5 files and does not create output directories.
#'
#' @param path Character scalar. POD5 file, directory, run tree, or manifest
#'   JSON file to plan from.
#' @param strategy One of `file-count`, `sample-label`, `elapsed-time`, or
#'   `read-count`.
#' @param files_per_chunk Positive integer used by `strategy = "file-count"`.
#' @param seconds_per_chunk Optional positive integer used by
#'   `strategy = "elapsed-time"`.
#' @param reads_per_chunk Optional positive integer used by
#'   `strategy = "read-count"`.
#'
#' @return A data frame with one row per planned chunk and the columns
#'   `schema_version`, `source`, `strategy`, `target`, `chunk_index`,
#'   `chunk_label`, `file_count`, `total_bytes`, `read_count`, `relative_paths`,
#'   and `warnings`.
#'
#' @examples
#' signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
#' run_dir <- file.path(tempdir(), "flounder-subdivide-plan-example")
#' unlink(run_dir, recursive = TRUE)
#' dir.create(run_dir, recursive = TRUE)
#' writeBin(c(signature, as.raw(rep(0, 16)), signature), file.path(run_dir, "a.pod5"))
#' pod5_subdivide_plan(run_dir, files_per_chunk = 1)
#'
#' @export
pod5_subdivide_plan <- function(
  path,
  strategy = c("file-count", "sample-label", "elapsed-time", "read-count"),
  files_per_chunk = 4L,
  seconds_per_chunk = NULL,
  reads_per_chunk = NULL
) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    .flounder_pod5_error("`path` must be a non-missing character scalar.",
      category = "path"
    )
  }
  strategy <- match.arg(strategy)
  files_per_chunk <- .flounder_positive_integer_scalar(
    files_per_chunk,
    "files_per_chunk"
  )
  seconds_per_chunk <- .flounder_optional_positive_integer_scalar(
    seconds_per_chunk,
    "seconds_per_chunk"
  )
  reads_per_chunk <- .flounder_optional_positive_integer_scalar(
    reads_per_chunk,
    "reads_per_chunk"
  )

  response <- .flounder_pod5_subdivide_plan_call(
    path,
    strategy,
    as.character(files_per_chunk),
    .flounder_optional_integer_label(seconds_per_chunk),
    .flounder_optional_integer_label(reads_per_chunk)
  )
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_pod5_error(
      response$error %||% "POD5 subdivision planning failed in the Rust extension.",
      category = response$category %||% "unknown"
    )
  }

  data <- response$data
  names(data) <- c(
    "schema_version",
    "source",
    "strategy",
    "target",
    "chunk_index",
    "chunk_label",
    "file_count",
    "total_bytes",
    "read_count",
    "relative_paths",
    "warnings"
  )
  data$schema_version <- as.integer(data$schema_version)
  data$chunk_index <- .flounder_nan_to_na(data$chunk_index)
  data$file_count <- as.integer(data$file_count)
  data$total_bytes <- as.numeric(data$total_bytes)
  data$read_count <- .flounder_nan_to_na(data$read_count)
  data$chunk_label <- .flounder_empty_to_na(data$chunk_label)
  data$relative_paths <- .flounder_empty_to_na(data$relative_paths)
  data$warnings <- .flounder_empty_to_na(data$warnings)
  data
}

.flounder_pod5_find_call <- function(path) {
  .Call("flounder_pod5_find", path, PACKAGE = "floundeR")
}

.flounder_pod5_verify_call <- function(path) {
  .Call("flounder_pod5_verify", path, PACKAGE = "floundeR")
}

.flounder_pod5_file_info_call <- function(path) {
  .Call("flounder_pod5_file_info", path, PACKAGE = "floundeR")
}

.flounder_pod5_folder_info_call <- function(path) {
  .Call("flounder_pod5_folder_info", path, PACKAGE = "floundeR")
}

.flounder_pod5_manifest_call <- function(path) {
  .Call("flounder_pod5_manifest", path, PACKAGE = "floundeR")
}

.flounder_pod5_compare_call <- function(left, right) {
  .Call("flounder_pod5_compare", left, right, PACKAGE = "floundeR")
}

.flounder_pod5_subdivide_plan_call <- function(
  path,
  strategy,
  files_per_chunk,
  seconds_per_chunk,
  reads_per_chunk
) {
  .Call(
    "flounder_pod5_subdivide_plan",
    path,
    strategy,
    files_per_chunk,
    seconds_per_chunk,
    reads_per_chunk,
    PACKAGE = "floundeR"
  )
}

.flounder_pod5_error <- function(
  message,
  category = "unknown",
  parent = NULL,
  call = NULL
) {
  subclass <- switch(
    category,
    path = "floundeR_pod5_path_error",
    format = "floundeR_pod5_format_error",
    schema = "floundeR_pod5_schema_error",
    integrity = "floundeR_pod5_integrity_error",
    "floundeR_pod5_unknown_error"
  )
  condition <- structure(
    list(message = message, parent = parent, call = call),
    class = c(
      subclass,
      "floundeR_pod5_error",
      "floundeR_error",
      "error",
      "condition"
    )
  )
  stop(condition)
}

`%||%` <- function(left, right) {
  if (is.null(left) || length(left) == 0L || is.na(left[[1L]])) {
    return(right)
  }

  left
}

.flounder_empty_to_na <- function(x) {
  x[x == ""] <- NA_character_
  x
}

.flounder_nan_to_na <- function(x) {
  x <- as.numeric(x)
  x[is.nan(x)] <- NA_real_
  x
}

.flounder_positive_integer_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1L || x != as.integer(x)) {
    .flounder_pod5_error(
      paste0("`", name, "` must be a positive integer scalar."),
      category = "format"
    )
  }
  as.integer(x)
}

.flounder_optional_positive_integer_scalar <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }
  .flounder_positive_integer_scalar(x, name)
}

.flounder_optional_integer_label <- function(x) {
  if (is.null(x)) {
    ""
  } else {
    as.character(x)
  }
}
