test_that("bam_check_index reports missing sidecar evidence", {
  skip_if_no_flounder_rust("BAM index check")

  path <- tempfile("flounder-index-", fileext = ".bam")
  bam_fixture(path)

  index <- floundeR::bam_check_index(path)

  expect_s3_class(index, "floundeR_bam_check_index")
  expect_named(index, c("status", "index", "candidates", "error"))
  expect_flounder_schema(index$status, "bam_check_index_status")
  expect_flounder_schema(index$index, "bam_check_index_index")
  expect_flounder_schema(index$candidates, "bam_check_index_candidates")
  expect_flounder_schema(index$error, "bam_validate_error")
  expect_true(index$status$ok)
  expect_false(index$index$present)
  expect_false(index$index$usable)
  expect_equal(index$index$compatibility, "absent")
  expect_equal(index$error$code, character())

  required <- floundeR::bam_check_index(path, require = TRUE)
  expect_false(required$status$ok)
  expect_equal(required$error$code, "missing_index")
})

test_that("bam_check_map reports mapped-read evidence", {
  skip_if_no_flounder_rust("BAM mapping check")

  path <- tempfile("flounder-map-", fileext = ".bam")
  bam_fixture(path)

  mapping <- floundeR::bam_check_map(path, sample_records = 10L)

  expect_s3_class(mapping, "floundeR_bam_check_map")
  expect_named(mapping, c("status", "index", "summary", "references"))
  expect_flounder_schema(mapping$status, "bam_check_map_status")
  expect_flounder_schema(mapping$index, "bam_check_map_index")
  expect_flounder_schema(mapping$summary, "bam_check_map_summary")
  expect_flounder_schema(mapping$references, "bam_check_map_references")
  expect_true(mapping$status$ok)
  expect_equal(mapping$status$mapping_status, "mapped")
  expect_equal(mapping$status$has_mapped_reads, "true")
  expect_equal(mapping$status$evidence_source, "scan")
  expect_false(mapping$index$present)
  expect_equal(mapping$summary$records_examined, 2)
  expect_equal(mapping$summary$mapped_records_observed, 1)
  expect_equal(mapping$summary$unmapped_records_observed, 1)
  expect_equal(mapping$references$name, "chr1")
})

test_that("bam_check_sort reports declared and observed sort evidence", {
  skip_if_no_flounder_rust("BAM sort check")

  path <- tempfile("flounder-sort-", fileext = ".bam")
  bam_fixture(path)

  sorting <- floundeR::bam_check_sort(path, sample_records = 10L)

  expect_s3_class(sorting, "floundeR_bam_check_sort")
  expect_flounder_schema(sorting, "bam_check_sort")
  expect_true(sorting$ok)
  expect_equal(sorting$command, "check_sort")
  expect_equal(sorting$declared_so, "coordinate")
  expect_true(sorting$observed_order %in% c("coordinate", "indeterminate"))
  expect_equal(sorting$records_examined, 2)
  expect_true(is.na(sorting$first_violation_reason))
})

test_that("bam_check_tag reports aux-tag absence evidence", {
  skip_if_no_flounder_rust("BAM tag check")

  path <- tempfile("flounder-tag-", fileext = ".bam")
  bam_fixture(path)

  tag <- floundeR::bam_check_tag(path, "RG", full_scan = TRUE)

  expect_s3_class(tag, "floundeR_bam_check_tag")
  expect_flounder_schema(tag, "bam_check_tag")
  expect_true(tag$ok)
  expect_equal(tag$command, "check_tag")
  expect_equal(tag$tag, "RG")
  expect_equal(tag$mode, "full_scan")
  expect_equal(tag$result, "absent_in_full_scan")
  expect_false(tag$tag_found)
  expect_equal(tag$records_examined, 2)
  expect_equal(tag$records_with_tag, 0)
  expect_true(tag$full_file_scanned)
  expect_true(is.na(tag$error_code))
})

test_that("BAM report-card checks validate arguments", {
  skip_if_no_flounder_rust("BAM report-card checks")

  expect_error(
    floundeR::bam_check_index(NA_character_),
    class = "floundeR_bam_path_error"
  )
  expect_error(
    floundeR::bam_check_map("missing.bam", sample_records = -1L),
    class = "floundeR_bam_argument_error"
  )
  expect_error(
    floundeR::bam_check_sort("missing.bam", strict = NA),
    class = "floundeR_bam_argument_error"
  )
  expect_error(
    floundeR::bam_check_tag("missing.bam", NA_character_),
    class = "floundeR_bam_argument_error"
  )
  expect_error(
    floundeR::bam_check_tag("missing.bam", "R"),
    class = "floundeR_bam_argument_error"
  )
})
