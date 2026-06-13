test_that("qc_yield_over_time returns a versioned tidy contract", {
  yield <- qc_yield_over_time(
    fixture_path("sequencing_summary_dorado.tsv"),
    resolution_minutes = 1)

  expect_flounder_schema(yield, "qc_yield_over_time")
  expect_equal(unique(yield$schema_version), "flounder.qc_yield_over_time.v1")
  expect_equal(sum(yield$read_count), 3L)
  expect_equal(sum(yield$bases), 96)
  expect_equal(max(yield$cumulative_read_count[yield$passes_filtering]), 2L)
  expect_equal(max(yield$cumulative_bases[yield$passes_filtering]), 68)
})

test_that("qc_read_length_distribution returns binned read-length data", {
  length_bins <- qc_read_length_distribution(
    fixture_path("sequencing_summary_dorado.tsv"),
    bins = 2)

  expect_flounder_schema(length_bins, "qc_read_length_distribution")
  expect_equal(
    unique(length_bins$schema_version),
    "flounder.qc_read_length_distribution.v1")
  expect_equal(sum(length_bins$read_count), 3L)
  expect_equal(sum(length_bins$bases), 96)
  expect_true(all(length_bins$read_length_bin_end >
    length_bins$read_length_bin_start))
})

test_that("qc_quality_distribution returns binned Q-score data", {
  quality_bins <- qc_quality_distribution(
    fixture_path("sequencing_summary_dorado.tsv"),
    bins = 2)

  expect_flounder_schema(quality_bins, "qc_quality_distribution")
  expect_equal(
    unique(quality_bins$schema_version),
    "flounder.qc_quality_distribution.v1")
  expect_equal(sum(quality_bins$read_count), 3L)
  expect_equal(sum(quality_bins$bases), 96)
  expect_true(all(quality_bins$qscore_bin_end >
    quality_bins$qscore_bin_start))
})

test_that("qc_channel_density returns per-channel yield and pass state", {
  density <- qc_channel_density(fixture_path("sequencing_summary_dorado.tsv"))

  expect_flounder_schema(density, "qc_channel_density")
  expect_equal(unique(density$schema_version), "flounder.qc_channel_density.v1")
  expect_equal(sum(density$read_count), 3L)
  expect_equal(sum(density$bases), 96)
  expect_equal(density$passed_read_count[density$channel == 11L], 1L)
  expect_equal(density$failed_read_count[density$channel == 13L], 1L)
})

test_that("qc_barcode_composition returns barcode balance data", {
  composition <- qc_barcode_composition(
    fixture_path("sequencing_summary_dorado.tsv"))

  expect_flounder_schema(composition, "qc_barcode_composition")
  expect_equal(
    unique(composition$schema_version),
    "flounder.qc_barcode_composition.v1")
  expect_equal(sum(composition$read_count), 3L)
  expect_equal(sum(composition$bases), 96)
  expect_equal(sum(composition$read_fraction), 1)
  expect_equal(sum(composition$bases_fraction), 1)
  expect_true("unclassified" %in% composition$barcode_arrangement)
})

test_that("qc_barcode_composition supports unbarcoded summaries", {
  composition <- qc_barcode_composition(flnDr("sequencing_summary.txt.bz2"))

  expect_flounder_schema(composition, "qc_barcode_composition")
  expect_equal(composition$barcode_arrangement, "unclassified")
  expect_equal(composition$read_fraction, 1)
})
