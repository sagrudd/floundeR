test_that("real-data QC evidence script dry-runs without network or basecalling", {
  script <- test_path("..", "..", "scripts", "build-real-data-qc-evidence.R")
  skip_if_not(file.exists(script))

  output_dir <- tempfile("flounder-real-data-qc-")
  env <- c(
    paste0("FLOUNDER_REAL_DATA_QC_DIR=", output_dir),
    "FLOUNDER_REAL_DATA_QC=false",
    "FLOUNDER_RUN_NETWORK_TESTS=false",
    "FLOUNDER_DGX_HOST=test-dgx-host"
  )
  result <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(script, "--dry-run"),
    env = env,
    stdout = TRUE,
    stderr = TRUE
  )

  expect_match(paste(result, collapse = "\n"), "Real-data QC evidence workflow")
  expect_true(file.exists(file.path(output_dir, "workflow-manifest.json")))
  expect_true(file.exists(file.path(output_dir, "analysis-status.tsv")))
  expect_true(file.exists(file.path(output_dir, "metadata", "source-metadata.tsv")))
  expect_true(file.exists(file.path(
    output_dir, "mnematikon-handoff", "mnematikon-handoff.tsv")))

  status <- read.delim(file.path(output_dir, "analysis-status.tsv"))
  expect_equal(status$status[[1]], "dry-run")
  expect_match(status$detail[[1]], "floundeR does not perform basecalling")

  handoff <- read.delim(file.path(
    output_dir, "mnematikon-handoff", "mnematikon-handoff.tsv"))
  expect_true(any(handoff$value == "Mnematikon, not floundeR"))
  expect_true(any(grepl("never provides direct basecalling", handoff$value)))
  expect_true(any(handoff$value == "test-dgx-host"))
})
