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
  qc_yield_over_time = c(
    "schema_version",
    "bin_start_seconds",
    "bin_end_seconds",
    "bin_start_minutes",
    "bin_end_minutes",
    "passes_filtering",
    "read_count",
    "bases",
    "cumulative_read_count",
    "cumulative_bases"
  ),
  qc_read_length_distribution = c(
    "schema_version",
    "read_length_bin_start",
    "read_length_bin_end",
    "passes_filtering",
    "read_count",
    "bases"
  ),
  qc_quality_distribution = c(
    "schema_version",
    "qscore_bin_start",
    "qscore_bin_end",
    "passes_filtering",
    "read_count",
    "bases"
  ),
  qc_channel_density = c(
    "schema_version",
    "channel",
    "read_count",
    "bases",
    "passed_read_count",
    "failed_read_count",
    "pass_fraction"
  ),
  qc_barcode_composition = c(
    "schema_version",
    "barcode_arrangement",
    "read_count",
    "bases",
    "passed_read_count",
    "failed_read_count",
    "read_fraction",
    "bases_fraction"
  ),
  qc_report_card = c(
    "schema_version",
    "check_id",
    "check_label",
    "status",
    "observed_value",
    "warn_threshold",
    "fail_threshold",
    "comparator",
    "details"
  ),
  ont_pod5_manifest = c(
    "file_name",
    "source_bucket",
    "source_key",
    "bytes",
    "last_modified_utc",
    "read_count",
    "state",
    "sha256"
  ),
  pod5_verify = c(
    "path",
    "size_bytes",
    "overall_status",
    "check",
    "category",
    "status",
    "detail"
  ),
  pod5_file_info = c(
    "path",
    "size_bytes",
    "flow_cell_id",
    "sequencing_kit",
    "read_count",
    "acquisition_start_utc",
    "duration_seconds",
    "pod5_version",
    "integrity_status",
    "integrity_reason"
  ),
  pod5_folder_info = c(
    "path",
    "pod5_file_count",
    "total_bytes",
    "total_reads",
    "flow_cell_ids",
    "sequencing_kits",
    "acquisition_start_utc",
    "acquisition_end_utc",
    "integrity_status",
    "integrity_reason",
    "failed_file_count",
    "verification_failed_count",
    "duplicate_file_names",
    "warnings"
  ),
  pod5_manifest = c(
    "schema_version",
    "source",
    "relative_path",
    "path",
    "size_bytes",
    "verification_status",
    "verification_failed_checks"
  ),
  pod5_compare = c(
    "status",
    "kind",
    "relative_path",
    "left_size_bytes",
    "right_size_bytes",
    "left_verification_status",
    "right_verification_status"
  ),
  pod5_subdivide_plan = c(
    "schema_version",
    "source",
    "strategy",
    "target",
    "chunk_index",
    "chunk_label",
    "file_count",
    "total_bytes",
    "read_count",
    "relative_paths",
    "warnings"
  )
)

expect_flounder_schema <- function(object, schema) {
  testthat::expect_true(
    schema %in% names(flounder_qc_schemas),
    info = paste("Unknown floundeR schema:", schema)
  )
  testthat::expect_named(object, flounder_qc_schemas[[schema]])
}
