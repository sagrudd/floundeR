test_that("minimal Rust capability function is callable through R wrappers", {
  capabilities <- floundeR::flounder_rust_capabilities()

  expect_type(capabilities, "list")
  expect_identical(capabilities$schema_version, "flounder.rust_capabilities.v1")
  expect_identical(capabilities$package, "floundeR")
  expect_identical(capabilities$rust_crate, "flounder-extendr")
  expect_identical(capabilities$rust_crate_version, "0.3.0")
  expect_true(capabilities$compiled_support)
  expect_identical(capabilities$pod5_tools, "linked")
  expect_identical(capabilities$bamana, "linked")
  expect_true(capabilities$porkchop %in% c("linked", "not_linked"))
  expect_identical(capabilities$grammateus, "not_linked")
})

test_that("Rust capability wrapper returns unavailable metadata when optional", {
  capabilities <- floundeR:::.flounder_rust_capabilities(
    required = FALSE,
    call = function() stop("native symbol missing", call. = FALSE)
  )

  expect_false(capabilities$compiled_support)
  expect_identical(capabilities$pod5_tools, "unavailable")
  expect_identical(capabilities$bamana, "unavailable")
  expect_identical(capabilities$porkchop, "unavailable")
  expect_identical(capabilities$grammateus, "unavailable")
  expect_match(capabilities$message, "native symbol missing")
})

test_that("Rust capability wrapper raises a typed unavailable condition", {
  expect_error(
    floundeR:::.flounder_rust_capabilities(
      required = TRUE,
      call = function() stop("native symbol missing", call. = FALSE)
    ),
    class = "floundeR_rust_unavailable"
  )
})
