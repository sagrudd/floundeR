test_that("payload-bearing Bamana failures remain QC evidence", {
  skip_if_no_flounder_rust("BAM envelope semantics")

  path <- tempfile("flounder-envelope-index-", fileext = ".bam")
  bam_fixture(path)

  index <- floundeR::bam_check_index(path, require = TRUE)
  expect_s3_class(index, "floundeR_bam_check_index")
  expect_false(index$status$ok)
  expect_equal(index$error$code, "missing_index")
  expect_match(index$error$message, "index", ignore.case = TRUE)

  truncated_path <- tempfile("flounder-envelope-eof-", fileext = ".bam")
  bam_fixture_without_eof(truncated_path)
  eof <- floundeR::bam_check_eof(truncated_path)
  expect_s3_class(eof, "floundeR_bam_check_eof")
  expect_false(eof$ok)
  expect_equal(eof$error_code, "truncated_file")
  expect_false(eof$complete)
})

test_that("no-payload Bamana failures become typed R conditions", {
  skip_if_no_flounder_rust("BAM error semantics")

  path_error <- tryCatch(
    floundeR::bam_summary("missing.bam"),
    error = identity
  )
  expect_s3_class(path_error, "floundeR_bam_path_error")
  expect_s3_class(path_error, "floundeR_bam_error")
  expect_named(path_error, c("message", "code", "detail", "hint", "parent", "call"))
  expect_equal(path_error$code, "file_not_found")

  tag_error <- tryCatch(
    floundeR::bam_check_tag("missing.bam", "R"),
    error = identity
  )
  expect_s3_class(tag_error, "floundeR_bam_argument_error")
  expect_s3_class(tag_error, "floundeR_bam_error")
  expect_named(tag_error, c("message", "code", "detail", "hint", "parent", "call"))
  expect_equal(tag_error$code, "invalid_tag")
  expect_match(tag_error$message, "tag", ignore.case = TRUE)
})
