pod5_signature_fixture <- function(path, middle_bytes = as.raw(rep(0, 16))) {
  signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
  writeBin(c(signature, middle_bytes, signature), path)
  path
}

test_that("pod5_folder_info summarises run-tree POD5 QC signals", {
  skip_if_no_flounder_rust("POD5 folder metadata")

  root <- tempfile("pod5-folder-info-")
  dir.create(file.path(root, "sample-a"), recursive = TRUE)
  dir.create(file.path(root, "sample-b"), recursive = TRUE)
  dir.create(file.path(root, "failed"), recursive = TRUE)
  pod5_signature_fixture(file.path(root, "sample-a", "a.pod5"))
  pod5_signature_fixture(file.path(root, "sample-b", "a.pod5"))
  writeBin(as.raw(c(0x8b, 0x50, 0x4f, 0x44)), file.path(root, "failed", "broken.pod5"))

  info <- floundeR::pod5_folder_info(root)

  expect_s3_class(info, "data.frame")
  expect_flounder_schema(info, "pod5_folder_info")
  expect_equal(nrow(info), 1L)
  expect_identical(info$path, root)
  expect_equal(info$pod5_file_count, 3L)
  expect_equal(info$total_bytes, 68)
  expect_true(is.na(info$total_reads))
  expect_identical(info$failed_file_count, 0L)
  expect_identical(info$verification_failed_count, 1L)
  expect_identical(info$duplicate_file_names, "a.pod5")
  expect_identical(info$integrity_status, "failed")
  expect_match(info$warnings, "duplicate POD5 file names detected")
  expect_match(info$warnings, "one or more files failed implemented verification checks")
})

test_that("pod5_manifest returns a versioned row-wise collection inventory", {
  skip_if_no_flounder_rust("POD5 manifests")

  root <- tempfile("pod5-manifest-")
  dir.create(file.path(root, "sample-a"), recursive = TRUE)
  dir.create(file.path(root, "failed"), recursive = TRUE)
  pod5_signature_fixture(file.path(root, "sample-a", "a.pod5"))
  writeBin(as.raw(c(0x8b, 0x50, 0x4f, 0x44)), file.path(root, "failed", "broken.pod5"))

  manifest <- floundeR::pod5_manifest(root)

  expect_s3_class(manifest, "data.frame")
  expect_flounder_schema(manifest, "pod5_manifest")
  expect_equal(nrow(manifest), 2L)
  expect_equal(unique(manifest$schema_version), 1L)
  expect_equal(unique(manifest$source), root)
  expect_setequal(manifest$relative_path, c(
    file.path("sample-a", "a.pod5"),
    file.path("failed", "broken.pod5")
  ))
  expect_true(all(manifest$size_bytes %in% c(4, 32)))
  expect_true(any(
    manifest$relative_path == file.path("sample-a", "a.pod5") &
      manifest$verification_status == "incomplete"
  ))
  expect_true(any(
    manifest$relative_path == file.path("failed", "broken.pod5") &
      manifest$verification_status == "failed" &
      manifest$verification_failed_checks > 0L
  ))
})

test_that("pod5_compare reports changed and missing manifest entries", {
  skip_if_no_flounder_rust("POD5 comparison")

  left <- tempfile("pod5-compare-left-")
  right <- tempfile("pod5-compare-right-")
  dir.create(left, recursive = TRUE)
  dir.create(right, recursive = TRUE)
  pod5_signature_fixture(file.path(left, "a.pod5"))
  pod5_signature_fixture(file.path(left, "b.pod5"))
  pod5_signature_fixture(file.path(right, "a.pod5"), as.raw(rep(1, 24)))
  pod5_signature_fixture(file.path(right, "c.pod5"))

  comparison <- floundeR::pod5_compare(left, right)

  expect_s3_class(comparison, "data.frame")
  expect_flounder_schema(comparison, "pod5_compare")
  expect_equal(unique(comparison$status), "different")
  expect_setequal(comparison$kind, c(
    "missing_from_right",
    "missing_from_left",
    "changed"
  ))
  expect_true(any(
    comparison$kind == "changed" &
      comparison$relative_path == "a.pod5" &
      comparison$left_size_bytes == 32 &
      comparison$right_size_bytes == 40
  ))
  expect_true(any(
    comparison$kind == "missing_from_right" &
      comparison$relative_path == "b.pod5"
  ))
  expect_true(any(
    comparison$kind == "missing_from_left" &
      comparison$relative_path == "c.pod5"
  ))
})

test_that("pod5_compare returns an explicit match row for equivalent inputs", {
  skip_if_no_flounder_rust("POD5 comparison")

  root <- tempfile("pod5-compare-match-")
  dir.create(root, recursive = TRUE)
  pod5_signature_fixture(file.path(root, "a.pod5"))

  comparison <- floundeR::pod5_compare(root, root)

  expect_flounder_schema(comparison, "pod5_compare")
  expect_equal(nrow(comparison), 1L)
  expect_identical(comparison$status, "match")
  expect_identical(comparison$kind, "match")
  expect_true(is.na(comparison$relative_path))
  expect_true(is.na(comparison$left_size_bytes))
  expect_true(is.na(comparison$right_size_bytes))
})

test_that("folder, manifest, and compare wrappers map path and format errors", {
  skip_if_no_flounder_rust("POD5 folder and manifest errors")

  missing <- file.path(tempdir(), "missing-flounder-pod5-folder")
  unlink(missing, recursive = TRUE)
  text_file <- tempfile("pod5-manifest-format-", fileext = ".txt")
  writeLines("not pod5", text_file)

  expect_error(
    floundeR::pod5_folder_info(missing),
    class = "floundeR_pod5_path_error"
  )
  expect_error(
    floundeR::pod5_manifest(text_file),
    class = "floundeR_pod5_format_error"
  )
  expect_error(
    floundeR::pod5_compare(missing, missing),
    class = "floundeR_pod5_path_error"
  )
})
