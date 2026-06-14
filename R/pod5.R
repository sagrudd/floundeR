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

.flounder_pod5_find_call <- function(path) {
  .Call("flounder_pod5_find", path, PACKAGE = "floundeR")
}

.flounder_pod5_verify_call <- function(path) {
  .Call("flounder_pod5_verify", path, PACKAGE = "floundeR")
}

.flounder_pod5_file_info_call <- function(path) {
  .Call("flounder_pod5_file_info", path, PACKAGE = "floundeR")
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
