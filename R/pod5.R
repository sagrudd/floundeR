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

.flounder_pod5_find_call <- function(path) {
  .Call("flounder_pod5_find", path, PACKAGE = "floundeR")
}

.flounder_pod5_error <- function(message, parent = NULL, call = NULL) {
  condition <- structure(
    list(message = message, parent = parent, call = call),
    class = c(
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
