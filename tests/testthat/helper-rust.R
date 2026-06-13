flounder_rust_available <- function() {
  identical(tolower(Sys.getenv("FLOUNDER_RUST_AVAILABLE", unset = "")), "true")
}

skip_if_no_flounder_rust <- function(feature = "compiled Rust support") {
  if (!flounder_rust_available()) {
    testthat::skip(
      paste(
        feature,
        "is unavailable.",
        "Set FLOUNDER_RUST_AVAILABLE=true when compiled Rust bindings are built."
      )
    )
  }
}
