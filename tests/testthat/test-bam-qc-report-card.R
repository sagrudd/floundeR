test_that("bam_qc_report_card returns BAM checks on the standard schema", {
  skip_if_no_flounder_rust("BAM QC report card")

  path <- tempfile("flounder-bam-card-", fileext = ".bam")
  bam_fixture(path)

  summary <- floundeR::bam_summary(path)
  index <- floundeR::bam_check_index(path)
  sorting <- floundeR::bam_check_sort(path, sample_records = 10L)
  tag <- floundeR::bam_check_tag(path, "RG", full_scan = TRUE)
  eof <- floundeR::bam_check_eof(path)
  validation <- floundeR::bam_validate(path)

  card <- floundeR::bam_qc_report_card(
    summary = summary,
    index = index,
    sorting = sorting,
    tags = tag,
    eof = eof,
    validation = validation,
    provenance_anomalies = 0,
    expected_tags = character(),
    thresholds = list(
      mapping_fraction_min = c(warn = 0.40, fail = 0.20),
      duplicate_fraction_max = c(warn = 0.10, fail = 0.20),
      qc_fail_fraction_max = c(warn = 0.05, fail = 0.10),
      mapq_zero_fraction_max = c(warn = 0.75, fail = 0.90),
      missing_or_unusable_index_max = c(warn = 0, fail = 0),
      stale_index_max = c(warn = 0, fail = 0),
      sorting_mismatch_max = c(warn = 0, fail = 0),
      missing_expected_tags_max = c(warn = 0, fail = 0),
      eof_absence_max = c(warn = 0, fail = 0),
      validation_finding_count_max = c(warn = 0, fail = 0),
      provenance_anomaly_count_max = c(warn = 0, fail = 0)
    ))

  expect_flounder_schema(card, "qc_report_card")
  expect_equal(unique(card$schema_version), "flounder.bam_qc_report_card.v1")
  expect_equal(nrow(card), 11L)
  expect_true(all(card$status %in% c("pass", "warn", "fail")))
  expect_equal(
    card$status[card$check_id == "bam_mapping_fraction"],
    "pass")
  expect_equal(
    card$observed_value[card$check_id == "bam_mapping_fraction"],
    0.5)
  expect_equal(
    card$status[card$check_id == "bam_missing_or_unusable_index"],
    "fail")
  expect_equal(
    card$status[card$check_id == "bam_eof_absence"],
    "pass")
  expect_equal(
    card$status[card$check_id == "bam_validation_findings"],
    "fail")
  expect_gte(
    card$observed_value[card$check_id == "bam_validation_findings"],
    1)
})

test_that("bam_qc_report_card reports missing tags and absent evidence", {
  skip_if_no_flounder_rust("BAM QC report card")

  path <- tempfile("flounder-bam-card-tags-", fileext = ".bam")
  bam_fixture(path)
  tag <- floundeR::bam_check_tag(path, "RG", full_scan = TRUE)

  card <- floundeR::bam_qc_report_card(
    tags = list(tag),
    expected_tags = c("RG", "MM"),
    provenance_anomalies = NULL)

  expect_equal(
    card$status[card$check_id == "bam_missing_expected_tags"],
    "fail")
  expect_equal(
    card$observed_value[card$check_id == "bam_missing_expected_tags"],
    2)
  expect_equal(
    card$status[card$check_id == "bam_mapping_fraction"],
    "warn")
  expect_true(is.na(
    card$observed_value[card$check_id == "bam_provenance_anomalies"]))
})

test_that("bam_qc_report_card validates expected tag inputs", {
  expect_error(
    floundeR::bam_qc_report_card(expected_tags = c("RG", NA_character_)),
    "`expected_tags` must be a character vector without NA values"
  )
})
