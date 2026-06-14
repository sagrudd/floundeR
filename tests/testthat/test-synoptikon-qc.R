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

test_that("write_synoptikon_qc output validates against the v1 JSON schema", {
  summary_path <- fixture_path("sequencing_summary_dorado.tsv")
  summary <- qc_run_summary(summary_path, source_id = "fixture-run")
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
  output <- tempfile(fileext = ".json")

  write_synoptikon_qc(
    output,
    run_summary = summary,
    report_cards = card,
    run_identity = list(
      run_id = "fixture-run",
      sample_id = "fixture-sample",
      flow_cell_id = "fixture-flow-cell"),
    governance_context = list(
      tenant_id = "fixture-tenant",
      governance_domain_id = "fixture-domain",
      mneion_work_request_id = "fixture-request"),
    input_provenance = list(
      input_id = "fixture-sequencing-summary",
      kind = "sequencing_summary",
      uri = "tests/testthat/fixtures/sequencing_summary_dorado.tsv",
      size_bytes = unname(file.info(summary_path)$size),
      sha256 = paste0(rep("0", 64), collapse = ""),
      last_modified_utc = "2026-06-14T06:21:36Z",
      metadata = list(source_id = "fixture-run")),
    limitations = list(
      limitation_id = "fixture-no-pod5",
      severity = "info",
      section = "pod5",
      message = "Fixture payload validates handoff shape without POD5 bytes."),
    payload_id = "fixture-schema-valid",
    generated_at_utc = "2026-06-14T06:21:36Z")

  payload <- jsonlite::fromJSON(output, simplifyVector = FALSE)
  expect_synoptikon_schema_valid(payload)
})

test_that("Synoptikon schema validation rejects incompatible payloads", {
  payload <- as_synoptikon_qc(
    payload_id = "fixture-schema-invalid",
    generated_at_utc = "2026-06-14T06:21:36Z")
  payload$schema_version <- "flounder.synoptikon_qc_payload.v2"

  schema <- jsonlite::fromJSON(
    synoptikon_payload_schema_path(),
    simplifyVector = FALSE)
  errors <- synoptikon_validate_schema_node(payload, schema, schema, "$")

  expect_true(any(grepl("schema_version expected const", errors, fixed = TRUE)))
})

test_that("local Synoptikon schema stays aligned with mnemosyne when present", {
  mnemosyne_schema <- synoptikon_mnemosyne_schema_path()
  if (is.null(mnemosyne_schema)) {
    skip(paste(
      "No canonical floundeR Synoptikon QC schema found in ../mnemosyne;",
      "checked:",
      paste(synoptikon_mnemosyne_schema_candidates(), collapse = ", ")
    ))
  }

  local_schema <- jsonlite::fromJSON(
    synoptikon_payload_schema_path(),
    simplifyVector = FALSE)
  expected_schema <- jsonlite::fromJSON(
    mnemosyne_schema,
    simplifyVector = FALSE)

  expect_equal(local_schema, expected_schema)
})
