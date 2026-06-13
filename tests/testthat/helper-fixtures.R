fixture_path <- function(...) {
  testthat::test_path("fixtures", ...)
}

fixture_lines <- function(...) {
  readLines(fixture_path(...), warn = FALSE)
}
