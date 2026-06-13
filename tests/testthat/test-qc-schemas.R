test_that("sequencing summary fixture matches QC schema", {
  sequencing_summary <- read.delim(
    fixture_path("sequencing_summary_dorado.tsv"),
    stringsAsFactors = FALSE
  )

  expect_flounder_schema(sequencing_summary, "sequencing_summary")
  expect_type(sequencing_summary$read_id, "character")
  expect_type(sequencing_summary$channel, "integer")
  expect_type(sequencing_summary$passes_filtering, "logical")
  expect_true(all(sequencing_summary$sequence_length_template > 0))
})

test_that("barcode summary fixture matches QC schema", {
  barcode_summary <- read.delim(
    fixture_path("barcoding_summary.tsv"),
    stringsAsFactors = FALSE
  )

  expect_flounder_schema(barcode_summary, "barcode_summary")
  expect_type(barcode_summary$read_id, "character")
  expect_type(barcode_summary$barcode_score, "double")
  expect_true(all(barcode_summary$barcode_score >= 0))
})

test_that("POD5 manifest fixture matches QC schema", {
  pod5_manifest <- read.delim(
    fixture_path("pod5_manifest.tsv"),
    stringsAsFactors = FALSE
  )

  expect_flounder_schema(pod5_manifest, "pod5_manifest")
  expect_equal(pod5_manifest$source_bucket, rep("ont-open-data", 2))
  expect_true(all(grepl("\\.pod5$", pod5_manifest$file_name)))
  expect_true(all(nchar(pod5_manifest$sha256) == 64))
})
