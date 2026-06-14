test_that("library_preparation_report_card returns requested checks", {
  kit_candidates <- data.frame(
    kit_id = c("SQK-LSK114", "SQK-RBK114"),
    normalized_score = c(0.92, 0.40),
    support_level = c("supported", "supported"),
    lifecycle_status = c("current", "current"),
    validation_status = c("validated", "validated")
  )
  adapter_primer <- data.frame(
    read_id = c("read-1", "read-2", "read-3"),
    motif_name = c("LSK_adapter", "primer_a", "LSK_adapter"),
    motif_kind = c("adapter", "primer", "adapter"),
    motif_family = c("ligation", "ligation", "ligation")
  )
  barcode <- data.frame(
    read_id = c("read-1", "read-2", "read-2", "read-3"),
    motif_name = c("barcode01", "barcode01", "barcode02", "barcode03"),
    motif_kind = rep("barcode", 4),
    motif_family = rep("barcode", 4)
  )
  cdna <- data.frame(
    read_id = c("read-1", "read-2", "read-3"),
    class = c("full_length", "partial", "unclassified"),
    classified = c(TRUE, TRUE, FALSE),
    full_length = c(TRUE, FALSE, FALSE)
  )

  card <- library_preparation_report_card(
    kit_candidates = kit_candidates,
    adapter_primer = adapter_primer,
    barcode = barcode,
    cdna = cdna,
    expected_kit_id = "SQK-LSK114",
    thresholds = list(
      unexpected_kit_max = c(warn = 0, fail = 0),
      adapter_burden_max = c(warn = 0.50, fail = 0.80),
      barcode_ambiguity_fraction_max = c(warn = 0.20, fail = 0.50),
      cdna_incomplete_fraction_max = c(warn = 0.50, fail = 0.80),
      kit_support_risk_max = c(warn = 0, fail = 1)
    ))

  expect_flounder_schema(card, "qc_report_card")
  expect_equal(
    unique(card$schema_version),
    "flounder.library_preparation_report_card.v1")
  expect_equal(nrow(card), 5L)
  expect_equal(
    card$status[card$check_id == "library_unexpected_kit"],
    "pass")
  expect_equal(
    card$status[card$check_id == "library_adapter_burden"],
    "warn")
  expect_equal(
    card$status[card$check_id == "library_barcode_ambiguity"],
    "warn")
  expect_equal(
    card$status[card$check_id == "library_cdna_incomplete_fraction"],
    "warn")
  expect_equal(
    card$status[card$check_id == "library_kit_support_risk"],
    "pass")
})

test_that("library_preparation_report_card flags unexpected and unsupported kits", {
  kit_candidates <- data.frame(
    kit_id = "SQK-RBK004",
    normalized_score = 0.88,
    support_level = "experimental",
    lifecycle_status = "retired",
    validation_status = "unknown"
  )

  card <- library_preparation_report_card(
    kit_candidates = kit_candidates,
    expected_kit_id = "SQK-LSK114")

  expect_equal(
    card$status[card$check_id == "library_unexpected_kit"],
    "fail")
  expect_equal(
    card$observed_value[card$check_id == "library_unexpected_kit"],
    1)
  expect_equal(
    card$status[card$check_id == "library_kit_support_risk"],
    "fail")
  expect_equal(
    card$observed_value[card$check_id == "library_kit_support_risk"],
    2)
})

test_that("library_preparation_report_card handles absent evidence as warnings", {
  card <- library_preparation_report_card()

  expect_equal(
    card$status[card$check_id == "library_unexpected_kit"],
    "pass")
  expect_equal(
    card$status[card$check_id == "library_adapter_burden"],
    "warn")
  expect_true(is.na(
    card$observed_value[card$check_id == "library_adapter_burden"]))
  expect_equal(
    card$status[card$check_id == "library_kit_support_risk"],
    "warn")
})

test_that("library_preparation_report_card validates inputs", {
  expect_error(
    library_preparation_report_card(expected_kit_id = c("a", "b")),
    "`expected_kit_id` must be a non-empty character scalar")
  expect_error(
    library_preparation_report_card(adapter_primer = list()),
    "`adapter_primer` must be a data frame or NULL")
})
