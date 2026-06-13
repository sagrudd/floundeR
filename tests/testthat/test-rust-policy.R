test_that("compiled Rust tests are disabled by default", {
  old <- Sys.getenv("FLOUNDER_RUST_AVAILABLE", unset = NA_character_)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("FLOUNDER_RUST_AVAILABLE")
    } else {
      Sys.setenv(FLOUNDER_RUST_AVAILABLE = old)
    }
  })

  Sys.unsetenv("FLOUNDER_RUST_AVAILABLE")
  expect_false(flounder_rust_available())

  Sys.setenv(FLOUNDER_RUST_AVAILABLE = "false")
  expect_false(flounder_rust_available())

  Sys.setenv(FLOUNDER_RUST_AVAILABLE = "true")
  expect_true(flounder_rust_available())
})
