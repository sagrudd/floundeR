#' @useDynLib floundeR, .registration = TRUE
#' @keywords internal
"_PACKAGE"

.flounder_rust_capabilities <- function() {
  payload <- .Call("flounder_rust_capabilities", PACKAGE = "floundeR")
  parts <- strsplit(payload, "|", fixed = TRUE)[[1]]

  list(
    schema_version = parts[[1]],
    package = "floundeR",
    rust_crate = parts[[2]],
    rust_crate_version = parts[[3]],
    extendr_api_version = "0.9.0",
    compiled_support = TRUE,
    pod5_tools = "not_linked",
    bamana = "not_linked",
    porkchop = "not_linked",
    grammateus = "not_linked"
  )
}
