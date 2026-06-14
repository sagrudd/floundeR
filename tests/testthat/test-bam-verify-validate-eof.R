test_that("bam_verify reports shallow BAM verification evidence", {
  skip_if_no_flounder_rust("BAM verification")

  path <- tempfile("flounder-verify-", fileext = ".bam")
  bam_fixture(path)

  verified <- floundeR::bam_verify(path)

  expect_s3_class(verified, "floundeR_bam_verify")
  expect_flounder_schema(verified, "bam_verify")
  expect_equal(verified$command, "verify")
  expect_true(verified$ok)
  expect_equal(verified$detected_format, "BAM")
  expect_equal(verified$container, "BGZF")
  expect_true(verified$is_bam)
  expect_true(verified$shallow_verified)
  expect_false(verified$deep_validated)
  expect_match(verified$checks_performed, "parsed_bam_header")
  expect_true(is.na(verified$error_code))
})

test_that("bam_validate returns structural validation evidence", {
  skip_if_no_flounder_rust("BAM validation")

  path <- tempfile("flounder-validate-", fileext = ".bam")
  bam_fixture(path)

  validation <- floundeR::bam_validate(path)

  expect_s3_class(validation, "floundeR_bam_validate")
  expect_named(validation, c("status", "summary", "findings", "error"))
  expect_flounder_schema(validation$status, "bam_validate_status")
  expect_flounder_schema(validation$summary, "bam_validate_summary")
  expect_flounder_schema(validation$findings, "bam_validate_findings")
  expect_flounder_schema(validation$error, "bam_validate_error")
  expect_true(validation$status$ok)
  expect_true(validation$status$valid)
  expect_equal(validation$status$mode, "full")
  expect_true(validation$summary$header_valid)
  expect_equal(validation$summary$records_examined, 2)
  expect_true(validation$summary$full_file_examined)
  expect_equal(validation$summary$errors, 0)
  expect_equal(validation$error$code, character())

  header_only <- floundeR::bam_validate(path, header_only = TRUE)
  expect_equal(header_only$status$mode, "header_only")
  expect_equal(header_only$summary$records_examined, 0)
  expect_false(header_only$summary$full_file_examined)
})

test_that("bam_check_eof reports present and missing EOF marker evidence", {
  skip_if_no_flounder_rust("BAM EOF check")

  complete_path <- tempfile("flounder-eof-complete-", fileext = ".bam")
  bam_fixture(complete_path)
  complete <- floundeR::bam_check_eof(complete_path)

  expect_s3_class(complete, "floundeR_bam_check_eof")
  expect_flounder_schema(complete, "bam_check_eof")
  expect_true(complete$ok)
  expect_true(complete$bgzf_eof_present)
  expect_true(complete$complete)
  expect_true(is.na(complete$error_code))

  truncated_path <- tempfile("flounder-eof-missing-", fileext = ".bam")
  bam_fixture_without_eof(truncated_path)
  missing <- floundeR::bam_check_eof(truncated_path)

  expect_flounder_schema(missing, "bam_check_eof")
  expect_false(missing$ok)
  expect_false(missing$bgzf_eof_present)
  expect_false(missing$complete)
  expect_equal(missing$error_code, "truncated_file")
})

test_that("bam_verify, bam_validate, and bam_check_eof validate arguments", {
  skip_if_no_flounder_rust("BAM verification and validation")

  expect_error(
    floundeR::bam_verify(NA_character_),
    class = "floundeR_bam_path_error"
  )
  expect_error(
    floundeR::bam_validate("missing.bam", max_errors = 0L),
    class = "floundeR_bam_argument_error"
  )
  expect_error(
    floundeR::bam_validate("missing.bam"),
    class = "floundeR_bam_path_error"
  )
  expect_error(
    floundeR::bam_check_eof(NA_character_),
    class = "floundeR_bam_path_error"
  )
  expect_error(
    floundeR::bam_check_eof("missing.bam"),
    class = "floundeR_bam_path_error"
  )
})
