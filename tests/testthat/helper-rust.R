flounder_rust_available <- function() {
  floundeR::flounder_rust_available()
}

skip_if_no_flounder_rust <- function(feature = "compiled Rust support") {
  floundeR::skip_if_no_flounder_rust(feature = feature)
}
