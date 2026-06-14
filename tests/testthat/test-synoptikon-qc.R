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

test_that("as_synoptikon_qc includes all required evidence sections", {
  summary <- qc_run_summary(
    fixture_path("sequencing_summary_dorado.tsv"),
    source_id = "fixture-run")
  flowcell <- qc_channel_density(
    fixture_path("sequencing_summary_dorado.tsv"))
  barcode <- qc_barcode_composition(
    fixture_path("sequencing_summary_dorado.tsv"))
  pod5 <- data.frame(
    schema_version = "flounder.pod5_integrity.v1",
    status = "warn",
    pod5_file_count = 2L,
    failed_file_count = 1L,
    stringsAsFactors = FALSE)
  bam <- data.frame(
    schema_version = "flounder.bam_qc_report_card.v1",
    status = "fail",
    check_id = "missing_index",
    stringsAsFactors = FALSE)
  library_preparation <- data.frame(
    schema_version = "flounder.library_preparation_qc.v1",
    status = "not_checked",
    evidence_kind = "porkchop_reserved",
    stringsAsFactors = FALSE)

  payload <- as_synoptikon_qc(
    run_summary = summary,
    flowcell = list(channel_density = flowcell),
    barcode = list(composition = barcode),
    pod5 = list(integrity = pod5),
    bam = list(report_card = bam),
    library_preparation = list(evidence = library_preparation),
    payload_id = "fixture-sections",
    generated_at_utc = "2026-06-14T06:03:36Z")

  sections <- payload$qc_sections
  expect_named(
    sections,
    c(
      "sequencing_summary", "flowcell", "barcode", "pod5", "bam",
      "library_preparation"))
  expect_true(all(vapply(sections, `[[`, logical(1), "present")))
  expect_equal(sections$flowcell$schema_version, "flounder.qc_channel_density.v1")
  expect_equal(
    sections$barcode$schema_version,
    "flounder.qc_barcode_composition.v1")
  expect_equal(sections$pod5$status, "warn")
  expect_equal(sections$bam$status, "fail")
  expect_equal(sections$library_preparation$status, "not_checked")
  expect_equal(sections$pod5$tables$integrity[[1]]$failed_file_count, 1L)
  expect_equal(sections$bam$tables$report_card[[1]]$check_id, "missing_index")
  expect_equal(
    sections$library_preparation$tables$evidence[[1]]$evidence_kind,
    "porkchop_reserved")
})
