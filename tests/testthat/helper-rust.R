flounder_rust_available <- function() {
  floundeR::flounder_rust_available()
}

skip_if_no_flounder_rust <- function(feature = "compiled Rust support") {
  floundeR::skip_if_no_flounder_rust(feature = feature)
}

skip_if_no_flounder_porkchop <- function(feature = "Porkchop integration") {
  capabilities <- skip_if_no_flounder_rust(feature = feature)
  if (!identical(capabilities$porkchop, "linked")) {
    testthat::skip(paste(
      feature,
      "is unavailable.",
      "Reinstall with CARGO_FEATURE_ARGS=--features=porkchop-integration and a local ../porkchop checkout."
    ))
  }
  invisible(capabilities)
}
