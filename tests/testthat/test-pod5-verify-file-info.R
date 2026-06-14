pod5_signature_fixture <- function(path, middle_bytes = as.raw(rep(0, 16))) {
  signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
  writeBin(c(signature, middle_bytes, signature), path)
  path
}

test_that("pod5_verify returns structured checks for a signature-valid candidate", {
  skip_if_no_flounder_rust("POD5 verification")

  path <- tempfile("pod5-verify-valid-", fileext = ".pod5")
  pod5_signature_fixture(path)

  verified <- floundeR::pod5_verify(path)

  expect_s3_class(verified, "data.frame")
  expect_flounder_schema(verified, "pod5_verify")
  expect_equal(unique(verified$path), path)
  expect_equal(unique(verified$size_bytes), 32)
  expect_equal(unique(verified$overall_status), "incomplete")
  expect_true(any(
    verified$check == "leading_signature" &
      verified$status == "passed"
  ))
  expect_true(any(
    verified$check == "trailing_signature" &
      verified$status == "passed"
  ))
  expect_true(any(
    verified$check == "combined_file_layout" &
      verified$status == "not_checked"
  ))
})

test_that("pod5_file_info returns filesystem metadata and unavailable integrity", {
  skip_if_no_flounder_rust("POD5 file metadata")

  path <- tempfile("pod5-file-info-", fileext = ".pod5")
  pod5_signature_fixture(path)

  info <- floundeR::pod5_file_info(path)

  expect_s3_class(info, "data.frame")
  expect_flounder_schema(info, "pod5_file_info")
  expect_equal(nrow(info), 1L)
  expect_identical(info$path, path)
  expect_equal(info$size_bytes, 32)
  expect_true(is.na(info$flow_cell_id))
  expect_true(is.na(info$sequencing_kit))
  expect_true(is.na(info$read_count))
  expect_true(is.na(info$acquisition_start_utc))
  expect_true(is.na(info$duration_seconds))
  expect_true(is.na(info$pod5_version))
  expect_identical(info$integrity_status, "unavailable")
  expect_match(info$integrity_reason, "POD5 parser backend not configured")
})

test_that("pod5_verify reports non-POD5 extension as a structured failure", {
  skip_if_no_flounder_rust("POD5 verification")

  path <- tempfile("pod5-verify-nonpod5-", fileext = ".txt")
  pod5_signature_fixture(path)

  verified <- floundeR::pod5_verify(path)

  expect_equal(unique(verified$overall_status), "failed")
  expect_true(any(
    verified$check == "extension" &
      verified$category == "extension" &
      verified$status == "failed"
  ))
})

test_that("pod5_file_info maps non-POD5 extension to a format condition", {
  skip_if_no_flounder_rust("POD5 file metadata")

  path <- tempfile("pod5-file-info-nonpod5-", fileext = ".txt")
  pod5_signature_fixture(path)

  expect_error(
    floundeR::pod5_file_info(path),
    class = "floundeR_pod5_format_error"
  )
})

test_that("pod5_verify reports truncated files as signature failures", {
  skip_if_no_flounder_rust("POD5 verification")

  path <- tempfile("pod5-verify-truncated-", fileext = ".pod5")
  writeBin(as.raw(c(0x8b, 0x50, 0x4f, 0x44)), path)

  verified <- floundeR::pod5_verify(path)

  expect_equal(unique(verified$overall_status), "failed")
  expect_true(all(
    verified$status[verified$category == "signature"] == "failed"
  ))
})

test_that("pod5_verify and pod5_file_info map missing paths to path conditions", {
  skip_if_no_flounder_rust("POD5 verification")

  path <- file.path(tempdir(), "missing-flounder-pod5-file.pod5")
  unlink(path)

  expect_error(
    floundeR::pod5_verify(path),
    class = "floundeR_pod5_path_error"
  )
  expect_error(
    floundeR::pod5_file_info(path),
    class = "floundeR_pod5_path_error"
  )
})

test_that("POD5 error helper maps schema and integrity categories to conditions", {
  expect_error(
    floundeR:::.flounder_pod5_error("schema problem", category = "schema"),
    class = "floundeR_pod5_schema_error"
  )
  expect_error(
    floundeR:::.flounder_pod5_error("integrity problem", category = "integrity"),
    class = "floundeR_pod5_integrity_error"
  )
})
