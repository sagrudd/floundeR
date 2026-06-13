test_that("qc_report_card returns a versioned pass/warn/fail schema", {
  card <- qc_report_card(
    fixture_path("sequencing_summary_dorado.tsv"),
    thresholds = list(
      pass_fraction_min = c(warn = 0.50, fail = 0.25),
      mean_qscore_min = c(warn = 10, fail = 8),
      n50_read_length_min = c(warn = 30, fail = 20),
      total_bases_min = c(warn = 80, fail = 50),
      channel_count_min = c(warn = 3, fail = 2),
      unclassified_fraction_max = c(warn = 0.25, fail = 0.50),
      barcode_max_fraction = c(warn = 0.40, fail = 0.80)
    ))

  expect_flounder_schema(card, "qc_report_card")
  expect_equal(unique(card$schema_version), "flounder.qc_report_card.v1")
  expect_equal(nrow(card), 7L)
  expect_true(all(card$status %in% c("pass", "warn", "fail")))
  expect_equal(
    card$status[card$check_id == "pass_fraction"],
    "pass")
  expect_equal(
    card$status[card$check_id == "unclassified_fraction"],
    "warn")
  expect_equal(
    card$status[card$check_id == "barcode_max_fraction"],
    "pass")
})

test_that("qc_report_card reports threshold failures", {
  card <- qc_report_card(
    fixture_path("sequencing_summary_dorado.tsv"),
    thresholds = list(
      pass_fraction_min = c(warn = 0.95, fail = 0.90),
      mean_qscore_min = c(warn = 20, fail = 15),
      n50_read_length_min = c(warn = 100, fail = 50),
      total_bases_min = c(warn = 1000, fail = 500),
      channel_count_min = c(warn = 10, fail = 5),
      unclassified_fraction_max = c(warn = 0.10, fail = 0.20),
      barcode_max_fraction = c(warn = 0.20, fail = 0.30)
    ))

  expect_true(all(card$status == "fail"))
  expect_true(all(grepl("failure threshold", card$details)))
})

test_that("qc_report_card handles missing observed values as warnings", {
  summary <- qc_run_summary(
    fixture_path("sequencing_summary_dorado.tsv"),
    source_id = "fixture-run")
  summary$mean_qscore <- NA_real_

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

  mean_qscore <- card[card$check_id == "mean_qscore", ]
  barcode_balance <- card[card$check_id == "barcode_max_fraction", ]

  expect_equal(mean_qscore$status, "warn")
  expect_true(grepl("missing", mean_qscore$details))
  expect_equal(barcode_balance$status, "warn")
  expect_true(is.na(barcode_balance$observed_value))
})

test_that("qc_report_card validates threshold pairs", {
  thresholds <- qc_report_card_thresholds()
  thresholds$mean_qscore_min <- c(warn = 10)

  expect_error(
    qc_report_card(
      fixture_path("sequencing_summary_dorado.tsv"),
      thresholds = thresholds),
    "must contain named warn and fail values")
})
