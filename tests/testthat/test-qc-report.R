library(floundeR)

local_without_grammateus_runtime <- function(code) {
  old_env <- Sys.getenv("GRAMMATEUS_HOME", unset = NA_character_)
  old_option <- getOption("floundeR.grammateus_home", default = NULL)
  on.exit({
    if (is.na(old_env)) {
      Sys.unsetenv("GRAMMATEUS_HOME")
    } else {
      Sys.setenv(GRAMMATEUS_HOME = old_env)
    }
    options(floundeR.grammateus_home = old_option)
  }, add = TRUE)
  Sys.setenv(GRAMMATEUS_HOME = tempfile("missing-grammateus-env-"))
  options(floundeR.grammateus_home = tempfile("missing-grammateus-option-"))
  force(code)
}

write_test_svg <- function(path) {
  writeLines(
    c(
      "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"120\" height=\"80\">",
      "  <rect width=\"120\" height=\"80\" fill=\"#ffffff\"/>",
      "  <circle cx=\"60\" cy=\"40\" r=\"20\" fill=\"#237A57\"/>",
      "</svg>"
    ),
    path,
    useBytes = TRUE
  )
  path
}

test_that("qc_report writes contract and manifest without Grammateus runtime", {
  local_without_grammateus_runtime({
    output_dir <- tempfile("qc-report-")
    dir.create(output_dir)
    svg_path <- write_test_svg(file.path(output_dir, "yield.svg"))
    figure <- grammateus_figure_from_file(
      path = svg_path,
      figure_id = "figure_yield_over_time",
      caption = "Cumulative sequencing yield over time.",
      alt_text = "Line chart showing cumulative bases over elapsed time.",
      produced_at_utc = "2026-06-14T11:42:09Z"
    )
    elements <- grammateus_qc_report_elements(
      run_metadata = data.frame(run_id = "run-report"),
      qc_summary = data.frame(metric = "read_count", value = 3),
      methods = "Fixture metrics were summarised by floundeR.",
      produced_at_utc = "2026-06-14T11:42:09Z",
      run_id = "run-report"
    )

    report <- qc_report(
      elements = elements,
      figures = list(figure),
      output_dir = output_dir,
      output = c("html", "pdf"),
      render = "if_available",
      report_id = "report_fixture_qc",
      title = "Fixture QC report",
      produced_at_utc = "2026-06-14T11:42:09Z",
      run_id = "run-report"
    )

    expect_s3_class(report, "flounder_qc_report")
    expect_equal(report$schema_version, "flounder.qc_report.v1")
    expect_true(file.exists(report$contract_path))
    expect_true(file.exists(report$manifest_path))
    expect_match(report$provenance$source_hash, "^sha256:[0-9a-f]{64}$")
    expect_equal(report$manifest$schema_version,
                 "flounder.qc_report_manifest.v1")
    expect_equal(report$manifest$element_count, elements$element_count)
    expect_equal(report$manifest$figure_count, 1L)
    expect_equal(report$manifest$theme$brand, "Mnemosyne Biosciences")
    expect_true(all(vapply(
      report$outputs,
      function(output) identical(output$status, "runtime_unavailable"),
      logical(1)
    )))

    contract <- jsonlite::fromJSON(report$contract_path, simplifyVector = FALSE)
    manifest <- jsonlite::fromJSON(report$manifest_path, simplifyVector = FALSE)
    expect_equal(contract$schema_version, "flounder.qc_report_contract.v1")
    expect_equal(contract$report_id, "report_fixture_qc")
    expect_equal(length(contract$figures), 1L)
    expect_equal(manifest$contract$sha256, report$manifest$contract$sha256)
    expect_match(manifest$contract$sha256, "^sha256:[0-9a-f]{64}$")
  })
})

test_that("qc_report supports manifest-only render policy", {
  output_dir <- tempfile("qc-report-manifest-only-")
  dir.create(output_dir)
  element <- grammateus_report_element(
    element_id = "table_qc_summary",
    element_type = "table",
    title = "QC summary",
    caption = "Run-level sequencing summary metrics.",
    data = data.frame(metric = "read_count", value = 3),
    produced_at_utc = "2026-06-14T11:42:09Z"
  )

  report <- qc_report(
    elements = element,
    output_dir = output_dir,
    output = "html",
    render = "never",
    report_id = "report_manifest_only",
    produced_at_utc = "2026-06-14T11:42:09Z"
  )

  expect_equal(length(report$outputs), 1L)
  expect_equal(report$outputs[[1L]]$format, "html")
  expect_equal(report$outputs[[1L]]$status, "not_requested")
  expect_true(file.exists(report$manifest_path))
  expect_true(file.exists(report$contract_path))
})

test_that("qc_report validates identifiers, formats, and required runtime", {
  output_dir <- tempfile("qc-report-validation-")
  dir.create(output_dir)
  element <- grammateus_report_element(
    element_id = "table_qc_summary",
    element_type = "table",
    title = "QC summary",
    caption = "Run-level sequencing summary metrics.",
    data = data.frame(metric = "read_count", value = 3)
  )

  expect_error(
    qc_report(
      elements = element,
      output_dir = output_dir,
      report_id = "bad_report"
    ),
    "report_id must start"
  )

  expect_error(
    qc_report(
      elements = element,
      output_dir = output_dir,
      output = "docx"
    ),
    "Unsupported qc_report output"
  )

  local_without_grammateus_runtime({
    expect_error(
      qc_report(
        elements = element,
        output_dir = output_dir,
        render = "require"
      ),
      class = "flounder_grammateus_runtime_unavailable"
    )
  })
})
