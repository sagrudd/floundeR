flounder_qc_schemas <- list(
  sequencing_summary = c(
    "read_id",
    "channel",
    "start_time",
    "duration",
    "passes_filtering",
    "sequence_length_template",
    "mean_qscore_template",
    "barcode_arrangement"
  ),
  barcode_summary = c(
    "read_id",
    "barcode_arrangement",
    "barcode_score",
    "barcode_front_id",
    "barcode_rear_id"
  ),
  qc_run_summary = c(
    "schema_version",
    "source_id",
    "read_count",
    "passed_read_count",
    "failed_read_count",
    "pass_fraction",
    "total_bases",
    "passed_bases",
    "failed_bases",
    "mean_read_length",
    "median_read_length",
    "n50_read_length",
    "mean_qscore",
    "median_qscore",
    "channel_count",
    "first_read_start_time",
    "last_read_start_time",
    "run_duration_seconds",
    "barcode_count",
    "unclassified_read_count"
  ),
  pod5_manifest = c(
    "file_name",
    "source_bucket",
    "source_key",
    "bytes",
    "last_modified_utc",
    "read_count",
    "state",
    "sha256"
  )
)

expect_flounder_schema <- function(object, schema) {
  testthat::expect_true(
    schema %in% names(flounder_qc_schemas),
    info = paste("Unknown floundeR schema:", schema)
  )
  testthat::expect_named(object, flounder_qc_schemas[[schema]])
}
