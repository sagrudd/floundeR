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
