test_that("minimal Rust capability function is callable from R", {
  capabilities <- floundeR:::.flounder_rust_capabilities()

  expect_type(capabilities, "list")
  expect_identical(capabilities$schema_version, "flounder.rust_capabilities.v1")
  expect_identical(capabilities$package, "floundeR")
  expect_identical(capabilities$rust_crate, "flounder-extendr")
  expect_identical(capabilities$rust_crate_version, "0.1.6")
  expect_true(capabilities$compiled_support)
  expect_identical(capabilities$pod5_tools, "not_linked")
  expect_identical(capabilities$bamana, "not_linked")
  expect_identical(capabilities$porkchop, "not_linked")
  expect_identical(capabilities$grammateus, "not_linked")
})
