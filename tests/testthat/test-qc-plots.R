library(floundeR)

test_that("QC plot helpers build sequencing-summary plot specs", {
  summary_path <- fixture_path("sequencing_summary_dorado.tsv")
  yield <- qc_yield_over_time(summary_path, resolution_minutes = 1)
  quality <- qc_quality_distribution(summary_path, bins = 2)
  read_lengths <- qc_read_length_distribution(summary_path, bins = 2)
  density <- qc_channel_density(summary_path)
  barcodes <- qc_barcode_composition(summary_path)

  specs <- list(
    qc_plot_yield_over_time(yield, produced_at_utc = "2026-06-14T11:57:09Z"),
    qc_plot_quality_distribution(quality, produced_at_utc = "2026-06-14T11:57:09Z"),
    qc_plot_read_length_distribution(read_lengths, produced_at_utc = "2026-06-14T11:57:09Z"),
    qc_plot_flowcell_density(density, produced_at_utc = "2026-06-14T11:57:09Z"),
    qc_plot_barcode_balance(barcodes, produced_at_utc = "2026-06-14T11:57:09Z")
  )
  ids <- vapply(specs, function(spec) spec$plot_id, character(1))

  expect_equal(ids, c(
    "plot_yield_over_time",
    "plot_quality_distribution",
    "plot_read_length_distribution",
    "plot_flowcell_density",
    "plot_barcode_balance"
  ))
  expect_true(all(vapply(
    specs,
    inherits,
    logical(1),
    what = "flounder_grammateus_plot_spec"
  )))
  expect_equal(specs[[1]]$plot_type, "line")
  expect_equal(specs[[2]]$plot_type, "stacked_bar")
  expect_equal(specs[[2]]$bar_value_semantics, "counts")
  expect_equal(specs[[5]]$bar_value_semantics, "rates")
  expect_match(specs[[1]]$provenance$source_hash, "^sha256:[0-9a-f]{64}$")
})

test_that("POD5 integrity plot helper summarises status-like evidence", {
  manifest <- data.frame(
    schema_version = "flounder.pod5_manifest.v1",
    relative_path = c("a.pod5", "b.pod5", "c.pod5"),
    verification_status = c("passed", "failed", "failed")
  )

  spec <- qc_plot_pod5_integrity(
    manifest,
    produced_at_utc = "2026-06-14T11:57:09Z"
  )

  expect_s3_class(spec, "flounder_grammateus_plot_spec")
  expect_equal(spec$plot_id, "plot_pod5_integrity")
  expect_equal(spec$plot_type, "bar")
  statuses <- vapply(spec$data$records, `[[`, character(1), "status")
  counts <- vapply(spec$data$records, `[[`, numeric(1), "file_count")
  expect_equal(stats::setNames(counts, statuses)[["failed"]], 2)
})

test_that("BAM plot helpers build mapping, MAPQ, and flag plot specs", {
  summary <- list(
    counts = data.frame(
      mapped_records = 8,
      unmapped_records = 2,
      primary_records = 9,
      secondary_records = 1,
      supplementary_records = 0,
      duplicate_records = 1,
      qc_fail_records = 0,
      paired_records = 4,
      properly_paired_records = 3,
      read1_records = 2,
      read2_records = 2
    ),
    mapq_histogram = data.frame(mapq = c(0L, 20L, 60L),
                                read_count = c(2, 3, 5))
  )
  class(summary) <- c("floundeR_bam_summary", "list")

  mapping <- qc_plot_bam_mapping_summary(
    summary,
    produced_at_utc = "2026-06-14T11:57:09Z"
  )
  mapq <- qc_plot_bam_mapq_distribution(
    summary,
    produced_at_utc = "2026-06-14T11:57:09Z"
  )
  flags <- qc_plot_bam_flag_summary(
    summary,
    produced_at_utc = "2026-06-14T11:57:09Z"
  )

  expect_equal(mapping$plot_id, "plot_bam_mapping_summary")
  expect_equal(mapq$plot_id, "plot_bam_mapq_distribution")
  expect_equal(flags$plot_id, "plot_bam_flag_summary")
  expect_equal(mapping$bar_value_semantics, "counts")
  expect_equal(length(mapping$data$records), 5L)
  expect_equal(length(mapq$data$records), 3L)
  expect_equal(length(flags$data$records), 6L)
})

test_that("QC plot helpers validate input columns", {
  expect_error(
    qc_plot_yield_over_time(data.frame(bin_end_minutes = 1)),
    "missing column"
  )
  expect_error(
    qc_plot_pod5_integrity(data.frame(path = "a.pod5")),
    "status or integrity"
  )
  expect_error(
    qc_plot_bam_mapping_summary(data.frame(mapped_records = 1)),
    "missing column"
  )
})
