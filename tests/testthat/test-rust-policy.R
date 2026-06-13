test_that("compiled Rust availability is detected from the native binding", {
  expect_true(flounder_rust_available())
  expect_true(floundeR::flounder_rust_available())
})

test_that("compiled Rust skip helper skips from capability metadata", {
  capabilities <- floundeR:::.flounder_unavailable_rust_capabilities(
    simpleError("native symbol missing")
  )

  expect_condition(
    floundeR:::.flounder_skip_if_no_flounder_rust(
      feature = "POD5 discovery",
      capabilities = capabilities
    ),
    class = "skip"
  )
})
