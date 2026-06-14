#' @useDynLib floundeR, .registration = TRUE, .fixes = "C_"
#' @keywords internal
"_PACKAGE"

#' Query compiled Rust support
#'
#' `flounder_rust_capabilities()` reports whether the compiled Rust extension is
#' available and which curated Rust-backed engines have been linked into the
#' package. `flounder_rust_available()` is a predicate for package code and
#' downstream tests. `skip_if_no_flounder_rust()` skips a `testthat` test when a
#' Rust-backed feature is unavailable.
#'
#' @param required Logical scalar. If `TRUE`, missing compiled Rust support is
#'   reported as a typed `floundeR_rust_unavailable` condition. If `FALSE`, a
#'   capability list with `compiled_support = FALSE` is returned.
#' @param feature Short feature label used in skip messages.
#'
#' @return
#' `flounder_rust_capabilities()` returns a named list with schema version,
#' package/crate metadata, `compiled_support`, and per-engine link status.
#' `flounder_rust_available()` returns `TRUE` or `FALSE`.
#' `skip_if_no_flounder_rust()` returns the capability list invisibly when Rust
#' is available, otherwise calls `testthat::skip()`.
#'
#' @export
flounder_rust_capabilities <- function(required = TRUE) {
  .flounder_rust_capabilities(required = required)
}

#' @rdname flounder_rust_capabilities
#' @export
flounder_rust_available <- function() {
  isTRUE(flounder_rust_capabilities(required = FALSE)$compiled_support)
}

#' @rdname flounder_rust_capabilities
#' @export
skip_if_no_flounder_rust <- function(feature = "compiled Rust support") {
  .flounder_skip_if_no_flounder_rust(feature = feature)
}

.flounder_rust_capabilities <- function(
  required = TRUE,
  call = .flounder_rust_capabilities_call
) {
  tryCatch(
    .flounder_parse_rust_capabilities(call()),
    error = function(error) {
      if (required) {
        .flounder_rust_unavailable(
          paste(
            "Compiled Rust support is unavailable.",
            "Install floundeR from source with Cargo and rustc available."
          ),
          parent = error
        )
      }

      .flounder_unavailable_rust_capabilities(error)
    }
  )
}

.flounder_rust_capabilities_call <- function() {
  .Call("flounder_rust_capabilities", PACKAGE = "floundeR")
}

.flounder_parse_rust_capabilities <- function(payload) {
  if (!is.character(payload) || length(payload) != 1L || is.na(payload)) {
    stop("Rust capability payload must be a non-missing character scalar.")
  }

  parts <- strsplit(payload, "|", fixed = TRUE)[[1]]
  if (length(parts) < 3L || length(parts) > 6L) {
    stop("Rust capability payload has an unsupported shape.")
  }

  pod5_tools <- if (length(parts) >= 4L && identical(parts[[4]], "pod5-tools")) {
    "linked"
  } else {
    "not_linked"
  }
  bamana <- if (length(parts) >= 5L && identical(parts[[5]], "bamana")) {
    "linked"
  } else {
    "not_linked"
  }
  porkchop <- if (length(parts) >= 6L && identical(parts[[6]], "porkchop")) {
    "linked"
  } else {
    "not_linked"
  }

  list(
    schema_version = parts[[1]],
    package = "floundeR",
    rust_crate = parts[[2]],
    rust_crate_version = parts[[3]],
    extendr_api_version = "0.9.0",
    compiled_support = TRUE,
    pod5_tools = pod5_tools,
    bamana = bamana,
    porkchop = porkchop,
    grammateus = "not_linked"
  )
}

.flounder_unavailable_rust_capabilities <- function(error) {
  list(
    schema_version = "flounder.rust_capabilities.v1",
    package = "floundeR",
    rust_crate = NA_character_,
    rust_crate_version = NA_character_,
    extendr_api_version = NA_character_,
    compiled_support = FALSE,
    pod5_tools = "unavailable",
    bamana = "unavailable",
    porkchop = "unavailable",
    grammateus = "unavailable",
    message = conditionMessage(error)
  )
}

.flounder_rust_unavailable <- function(message, parent = NULL, call = NULL) {
  condition <- structure(
    list(message = message, parent = parent, call = call),
    class = c(
      "floundeR_rust_unavailable",
      "floundeR_error",
      "error",
      "condition"
    )
  )
  stop(condition)
}

.flounder_skip_if_no_flounder_rust <- function(
  feature = "compiled Rust support",
  capabilities = flounder_rust_capabilities(required = FALSE)
) {
  if (isTRUE(capabilities$compiled_support)) {
    return(invisible(capabilities))
  }

  message <- paste(
    feature,
    "is unavailable.",
    "Install floundeR from source with Cargo and rustc available."
  )

  if (!requireNamespace("testthat", quietly = TRUE)) {
    .flounder_rust_unavailable(message)
  }

  testthat::skip(message)
}
