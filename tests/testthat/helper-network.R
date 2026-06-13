flounder_network_tests_enabled <- function() {
  identical(tolower(Sys.getenv("FLOUNDER_RUN_NETWORK_TESTS", unset = "")), "true")
}

skip_if_no_flounder_network <- function() {
  testthat::skip_on_cran()
  if (!flounder_network_tests_enabled()) {
    testthat::skip(
      paste(
        "Network-dependent floundeR tests are opt-in.",
        "Set FLOUNDER_RUN_NETWORK_TESTS=true to run them."
      )
    )
  }
}
