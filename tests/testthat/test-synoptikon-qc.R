test_that("as_synoptikon_qc builds the v1 payload shape", {
  summary <- qc_run_summary(
    fixture_path("sequencing_summary_dorado.tsv"),
    source_id = "fixture-run")
  card <- qc_report_card(
    summary,
    thresholds = list(
      pass_fraction_min = c(warn = 0.50, fail = 0.25),
      mean_qscore_min = c(warn = 10, fail = 8),
      n50_read_length_min = c(warn = 30, fail = 20),
      total_bases_min = c(warn = 80, fail = 50),
      channel_count_min = c(warn = 3, fail = 2),
      unclassified_fraction_max = c(warn = 0.25, fail = 0.50),
      barcode_max_fraction = c(warn = 0.40, fail = 0.80)
    ))

  payload <- as_synoptikon_qc(
    run_summary = summary,
    report_cards = card,
    run_identity = list(run_id = "fixture-run"),
    payload_id = "fixture-payload",
    generated_at_utc = "2026-06-14T05:43:36Z")

  expect_named(payload, flounder_qc_schemas$synoptikon_qc_payload)
  expect_equal(payload$schema_version, "flounder.synoptikon_qc_payload.v1")
  expect_equal(payload$payload_id, "fixture-payload")
  expect_equal(payload$run_identity$run_id, "fixture-run")
  expect_true(payload$qc_sections$sequencing_summary$present)
  expect_equal(
    payload$qc_sections$sequencing_summary$schema_version,
    "flounder.qc_run_summary.v1")
  expect_equal(
    payload$qc_sections$sequencing_summary$summary$read_count,
    3L)
  expect_length(payload$qc_sections$sequencing_summary$tables$run_summary, 1)
  expect_length(payload$report_cards, 1)
  expect_equal(payload$report_cards[[1]]$scope, "run")
  expect_length(payload$report_cards[[1]]$checks, 7)
})

test_that("write_synoptikon_qc writes JSON", {
  summary <- qc_run_summary(
    fixture_path("sequencing_summary_dorado.tsv"),
    source_id = "fixture-run")
  output <- tempfile(fileext = ".json")

  path <- write_synoptikon_qc(
    output,
    run_summary = summary,
    payload_id = "fixture-payload",
    generated_at_utc = "2026-06-14T05:43:36Z")

  expect_identical(path, output)
  parsed <- jsonlite::fromJSON(output, simplifyVector = FALSE)
  expect_equal(parsed$schema_version, "flounder.synoptikon_qc_payload.v1")
  expect_equal(parsed$payload_id, "fixture-payload")
  expect_true(parsed$qc_sections$sequencing_summary$present)
  expect_length(parsed$qc_sections$flowcell$summary, 0)
})
