test_that("qc_run_summary returns the v1 run summary schema", {
  summary <- qc_run_summary(
    fixture_path("sequencing_summary_dorado.tsv"),
    source_id = "fixture-run")

  expect_flounder_schema(summary, "qc_run_summary")
  expect_equal(summary$schema_version, "flounder.qc_run_summary.v1")
  expect_equal(summary$source_id, "fixture-run")
  expect_equal(summary$read_count, 3L)
  expect_equal(summary$passed_read_count, 2L)
  expect_equal(summary$failed_read_count, 1L)
  expect_equal(summary$pass_fraction, 2 / 3)
  expect_equal(summary$total_bases, 96)
  expect_equal(summary$passed_bases, 68)
  expect_equal(summary$failed_bases, 28)
  expect_equal(summary$mean_read_length, 32)
  expect_equal(summary$median_read_length, 32)
  expect_equal(summary$n50_read_length, 32)
  expect_equal(summary$mean_qscore, 11.033333, tolerance = 0.000001)
  expect_equal(summary$median_qscore, 11.8)
  expect_equal(summary$channel_count, 3L)
  expect_equal(summary$first_read_start_time, 0.125)
  expect_equal(summary$last_read_start_time, 5.25)
  expect_equal(summary$run_duration_seconds, 9.25)
  expect_equal(summary$barcode_count, 2L)
  expect_equal(summary$unclassified_read_count, 1L)
})

test_that("qc_run_summary accepts SequencingSummary objects and tibbles", {
  seqsum <- SequencingSummary$new(fixture_path("sequencing_summary_dorado.tsv"))
  object_summary <- qc_run_summary(seqsum, source_id = "object")
  tibble_summary <- qc_run_summary(seqsum$as_tibble(), source_id = "tibble")

  expect_equal(object_summary$read_count, tibble_summary$read_count)
  expect_equal(object_summary$total_bases, tibble_summary$total_bases)
  expect_equal(object_summary$source_id, "object")
  expect_equal(tibble_summary$source_id, "tibble")
})

test_that("qc_run_summary supports unbarcoded sequencing summaries", {
  summary <- qc_run_summary(flnDr("sequencing_summary.txt.bz2"))

  expect_flounder_schema(summary, "qc_run_summary")
  expect_true(summary$read_count > 0)
  expect_equal(summary$barcode_count, 0L)
  expect_equal(summary$unclassified_read_count, 0L)
})

test_that("qc_run_summary reports missing contract columns", {
  seqsum <- SequencingSummary$new(
    fixture_path("sequencing_summary_dorado.tsv"))$as_tibble()
  seqsum$channel <- NULL

  expect_error(
    qc_run_summary(seqsum),
    "requires sequencing-summary column\\(s\\): channel")
})
