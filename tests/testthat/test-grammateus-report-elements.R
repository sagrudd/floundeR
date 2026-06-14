library(floundeR)

local_fixture_path <- function(...) {
  candidates <- c(
    file.path("fixtures", ...),
    file.path("tests", "testthat", "fixtures", ...)
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0L) {
    return(existing[[1L]])
  }
  candidates[[1L]]
}

test_that("grammateus_report_element builds table payloads with provenance", {
  element <- grammateus_report_element(
    element_id = "table_qc_summary",
    element_type = "table",
    title = "QC summary",
    caption = "Run-level sequencing summary metrics.",
    data = data.frame(metric = "read_count", value = 3),
    produced_at_utc = as.POSIXct("2026-06-14 11:12:09", tz = "UTC"),
    run_id = "run-elements-001"
  )

  expect_s3_class(element, "flounder_grammateus_report_element")
  expect_equal(element$schema_version,
               "flounder.grammateus_report_element.v1")
  expect_equal(element$element_id, "table_qc_summary")
  expect_equal(element$element_type, "table")
  expect_equal(element$payload$kind, "table")
  expect_equal(element$payload$row_count, 1L)
  expect_equal(element$payload$columns, c("metric", "value"))
  expect_match(element$payload$data_sha256, "^sha256:[0-9a-f]{64}$")
  expect_match(element$provenance$source_hash, "^sha256:[0-9a-f]{64}$")
  expect_equal(element$provenance$produced_at_utc, "2026-06-14T11:12:09Z")
  expect_equal(element$provenance$run_id, "run-elements-001")
})

test_that("grammateus_qc_report_elements covers required QC report sections", {
  summary_path <- local_fixture_path("sequencing_summary_dorado.tsv")
  summary <- qc_run_summary(summary_path, source_id = "run-elements")
  yield <- qc_yield_over_time(summary_path, resolution_minutes = 1)
  quality <- qc_quality_distribution(summary_path, bins = 3)
  flowcell <- qc_channel_density(summary_path)
  barcode <- qc_barcode_composition(summary_path)
  card <- qc_report_card(summary, barcode_composition = barcode)
  pod5 <- data.frame(schema_version = "flounder.pod5_integrity.v1",
                     status = "pass", pod5_file_count = 1L)
  bam_summary <- data.frame(schema_version = "flounder.bam_summary.v1",
                            mapped_fraction = 0.92)
  bam_validation <- data.frame(schema_version = "flounder.bam_validation.v1",
                               status = "pass")
  bam_index <- data.frame(schema_version = "flounder.bam_index.v1",
                          status = "present")
  bam_sort <- data.frame(schema_version = "flounder.bam_sort.v1",
                         status = "coordinate")
  library_prep <- data.frame(
    schema_version = "flounder.library_preparation_qc.v1",
    evidence_kind = "kit_candidate",
    status = "pass"
  )
  provenance <- data.frame(input_id = "seqsum", kind = "sequencing_summary")

  bundle <- grammateus_qc_report_elements(
    run_metadata = data.frame(run_id = "run-elements"),
    qc_summary = as.data.frame(summary),
    flowcell_density = flowcell,
    yield_over_time = yield,
    quality_distribution = quality,
    barcode_balance = barcode,
    pod5_integrity = pod5,
    bam_alignment_summary = bam_summary,
    bam_validation = bam_validation,
    bam_index = bam_index,
    bam_sort = bam_sort,
    library_preparation = library_prep,
    report_card_findings = card,
    methods = "Sequencing-summary fixture metrics were summarised by floundeR.",
    limitations = data.frame(section = "pod5", message = "Fixture POD5 only."),
    appendices = list(extra_metrics = data.frame(metric = "x", value = 1)),
    provenance = provenance,
    produced_at_utc = "2026-06-14T11:12:09Z",
    run_id = "run-elements"
  )

  expect_s3_class(bundle, "flounder_grammateus_report_element_bundle")
  expected <- c(
    "run_metadata", "qc_summary", "flowcell_density", "yield_over_time",
    "quality_distribution", "barcode_balance", "pod5_integrity",
    "bam_alignment_summary", "bam_validation", "bam_index", "bam_sort",
    "library_preparation", "report_card_findings", "methods", "limitations",
    "provenance", "appendix_extra_metrics"
  )
  expect_named(bundle$elements, expected)
  expect_equal(bundle$element_count, length(expected))
  expect_equal(bundle$elements$yield_over_time$element_id,
               "table_yield_over_time")
  expect_equal(bundle$elements$library_preparation$title,
               "Library-preparation evidence")
  expect_equal(bundle$elements$methods$element_type, "methods")
  expect_equal(bundle$elements$limitations$element_type, "limitations")
  expect_equal(bundle$elements$appendix_extra_metrics$element_type, "appendix")
  expect_true(all(vapply(
    bundle$elements,
    function(element) grepl("^sha256:[0-9a-f]{64}$",
                            element$provenance$source_hash),
    logical(1)
  )))
})

test_that("grammateus report elements validate identifiers and captions", {
  expect_error(
    grammateus_report_element(
      element_id = "bad_qc_summary",
      element_type = "table",
      title = "Bad table",
      caption = "Bad identifier.",
      data = data.frame(x = 1)
    ),
    "element_id must start"
  )

  expect_error(
    grammateus_report_element(
      element_id = "table_bad",
      element_type = "table",
      title = "Bad table",
      caption = "",
      data = data.frame(x = 1)
    ),
    "caption must be"
  )

  expect_error(
    grammateus_qc_report_elements(
      appendices = stats::setNames(list("bad appendix"), "")
    ),
    "appendices must be a named list"
  )
})
