test_that("network tests are disabled by default", {
  old <- Sys.getenv("FLOUNDER_RUN_NETWORK_TESTS", unset = NA_character_)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("FLOUNDER_RUN_NETWORK_TESTS")
    } else {
      Sys.setenv(FLOUNDER_RUN_NETWORK_TESTS = old)
    }
  })

  Sys.unsetenv("FLOUNDER_RUN_NETWORK_TESTS")
  expect_false(flounder_network_tests_enabled())

  Sys.setenv(FLOUNDER_RUN_NETWORK_TESTS = "false")
  expect_false(flounder_network_tests_enabled())

  Sys.setenv(FLOUNDER_RUN_NETWORK_TESTS = "true")
  expect_true(flounder_network_tests_enabled())
})
