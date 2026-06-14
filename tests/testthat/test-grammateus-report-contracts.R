library(floundeR)

local_without_grammateus_runtime_contracts <- function(code) {
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

write_contract_svg <- function(path) {
  writeLines(
    c(
      "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"120\" height=\"80\">",
      "  <rect width=\"120\" height=\"80\" fill=\"#ffffff\"/>",
      "  <path d=\"M10 60 L55 35 L110 20\" stroke=\"#237A57\"",
      "        stroke-width=\"6\" fill=\"none\"/>",
      "</svg>"
    ),
    path,
    useBytes = TRUE
  )
  path
}

expect_sha256 <- function(value, label = deparse(substitute(value))) {
  expect_true(is.character(value), label = label)
  expect_length(value, 1L)
  expect_match(value, "^sha256:[0-9a-f]{64}$")
}

expect_non_empty_scalar <- function(value, label = deparse(substitute(value))) {
  expect_true(is.character(value), label = label)
  expect_length(value, 1L)
  expect_true(nzchar(trimws(value)), label = label)
}

expect_report_element_contract <- function(element) {
  expect_s3_class(element, "flounder_grammateus_report_element")
  expect_equal(element$schema_version, "flounder.grammateus_report_element.v1")
  expect_match(element$element_id, "^[a-z0-9]+_[a-z0-9_]+$")
  prefix_ok <- startsWith(element$element_id, paste0(element$element_type, "_"))
  prefix_ok <- prefix_ok || (element$element_type == "table" &&
    startsWith(element$element_id, "table_"))
  expect_true(prefix_ok, label = element$element_id)
  expect_non_empty_scalar(element$title, paste(element$element_id, "title"))
  if (element$element_type %in%
      c("table", "methods", "limitations", "appendix", "provenance")) {
    expect_non_empty_scalar(
      element$caption,
      paste(element$element_id, "caption")
    )
  }
  expect_true(is.list(element$payload))
  expect_true(element$payload$kind %in% c("table", "list", "text"))
  if (element$payload$kind %in% c("table", "list")) {
    expect_sha256(
      element$payload$data_sha256,
      paste(element$element_id, "payload data hash")
    )
  }
  expect_sha256(
    element$provenance$source_hash,
    paste(element$element_id, "provenance hash")
  )
  expect_non_empty_scalar(element$provenance$produced_by)
  expect_non_empty_scalar(element$provenance$producer_version)
  expect_non_empty_scalar(element$provenance$produced_at_utc)
}

expect_figure_contract <- function(figure) {
  expect_s3_class(figure, "flounder_grammateus_figure")
  expect_equal(figure$schema_version, "flounder.grammateus_figure.v1")
  expect_match(figure$figure_id, "^figure_[a-z0-9_]+$")
  expect_non_empty_scalar(figure$caption)
  expect_non_empty_scalar(figure$alt_text)
  expect_equal(figure$source$kind, "image")
  expect_true(figure$source$format %in% c("png", "svg"))
  expect_true(file.exists(figure$source$path))
  expect_sha256(figure$source$checksum)
  expect_gt(figure$source$width_px, 0)
  expect_gt(figure$source$height_px, 0)
  expect_sha256(figure$provenance$source_hash)
}

test_that("standard QC report element bundles satisfy Grammateus contracts", {
  bundle <- grammateus_qc_report_elements(
    run_metadata = data.frame(run_id = "run-contract", sample_id = "sample-a"),
    qc_summary = data.frame(metric = "read_count", value = 3),
    yield_over_time = data.frame(minutes = c(0, 1), bases = c(10, 25)),
    pod5_integrity = data.frame(status = "pass", pod5_file_count = 1L),
    bam_alignment_summary = data.frame(category = "mapped", read_count = 2L),
    library_preparation = data.frame(
      evidence_kind = "kit_candidate",
      status = "pass",
      score_type = "heuristic_score"
    ),
    report_card_findings = data.frame(
      check_id = "read_count",
      status = "pass",
      details = "Fixture check."
    ),
    methods = "Fixture data were summarised with floundeR.",
    limitations = "Tiny fixture data are not representative of production runs.",
    provenance = data.frame(input_id = "seqsum", sha256 = paste0(
      "sha256:",
      paste(rep("a", 64), collapse = "")
    )),
    appendices = list(extra_metrics = data.frame(metric = "n50", value = 42)),
    produced_at_utc = "2026-06-14T12:12:09Z",
    run_id = "run-contract"
  )

  expect_s3_class(bundle, "flounder_grammateus_report_element_bundle")
  expect_equal(bundle$schema_version,
               "flounder.grammateus_report_element_bundle.v1")
  expect_equal(bundle$element_count, length(bundle$elements))
  element_ids <- vapply(
    bundle$elements,
    function(element) element$element_id,
    character(1)
  )
  expect_identical(anyDuplicated(element_ids), 0L)
  invisible(lapply(bundle$elements, expect_report_element_contract))
})

test_that("figures and plot specs carry required captions, alt text, and hashes", {
  output_dir <- tempfile("grammateus-contract-")
  dir.create(output_dir)
  svg_path <- write_contract_svg(file.path(output_dir, "yield.svg"))
  figure <- grammateus_figure_from_file(
    path = svg_path,
    figure_id = "figure_yield_over_time",
    caption = "Cumulative sequencing yield over time.",
    alt_text = "Line chart showing cumulative base yield over elapsed time.",
    methods_note = "Generated from sequencing summary fixture records.",
    source_hash = paste(rep("b", 64), collapse = ""),
    produced_at_utc = "2026-06-14T12:12:09Z",
    run_id = "run-contract"
  )
  expect_figure_contract(figure)
  expect_equal(figure$provenance$source_hash, paste0(
    "sha256:",
    paste(rep("b", 64), collapse = "")
  ))

  plot <- qc_plot_yield_over_time(
    data.frame(
      bin_end_minutes = c(1, 2),
      cumulative_bases = c(100, 250),
      passes_filtering = c(TRUE, TRUE)
    ),
    produced_at_utc = "2026-06-14T12:12:09Z",
    run_id = "run-contract"
  )
  expect_s3_class(plot, "flounder_grammateus_plot_spec")
  expect_equal(plot$schema_version, "flounder.grammateus_plot.v1")
  expect_match(plot$plot_id, "^plot_[a-z0-9_]+$")
  expect_non_empty_scalar(plot$caption)
  expect_equal(plot$data$kind, "inline_tidy")
  expect_length(plot$data$records, 2L)
  expect_sha256(plot$provenance$source_hash)
})

test_that("written qc_report contracts preserve stable identifiers and hashes", {
  local_without_grammateus_runtime_contracts({
    output_dir <- tempfile("qc-report-contract-")
    dir.create(output_dir)
    svg_path <- write_contract_svg(file.path(output_dir, "yield.svg"))
    figure <- grammateus_figure_from_file(
      path = svg_path,
      figure_id = "figure_yield_over_time",
      caption = "Cumulative sequencing yield over time.",
      alt_text = "Line chart showing cumulative base yield over elapsed time.",
      produced_at_utc = "2026-06-14T12:12:09Z",
      run_id = "run-contract"
    )
    elements <- grammateus_qc_report_elements(
      run_metadata = data.frame(run_id = "run-contract"),
      qc_summary = data.frame(metric = "read_count", value = 3),
      methods = "Fixture data were summarised with floundeR.",
      produced_at_utc = "2026-06-14T12:12:09Z",
      run_id = "run-contract"
    )

    report <- qc_report(
      elements = elements,
      figures = figure,
      output_dir = output_dir,
      output = c("html", "pdf"),
      render = "if_available",
      report_id = "report_contract_regression",
      title = "Contract regression QC report",
      produced_at_utc = "2026-06-14T12:12:09Z",
      run_id = "run-contract"
    )

    expect_sha256(report$provenance$source_hash)
    expect_sha256(report$manifest$contract$sha256)
    expect_true(file.exists(report$contract_path))
    expect_true(file.exists(report$manifest_path))
    expect_equal(
      report$manifest$contract$sha256,
      paste0("sha256:", unname(tools::sha256sum(report$contract_path)))
    )

    contract <- jsonlite::fromJSON(report$contract_path, simplifyVector = FALSE)
    manifest <- jsonlite::fromJSON(report$manifest_path, simplifyVector = FALSE)
    expect_equal(contract$schema_version, "flounder.qc_report_contract.v1")
    expect_equal(contract$report_id, "report_contract_regression")
    expect_equal(contract$themed_report$element_count, 3L)
    expect_equal(length(contract$figures), 1L)
    expect_equal(contract$figures[[1L]]$figure_id, "figure_yield_over_time")
    expect_non_empty_scalar(contract$figures[[1L]]$caption)
    expect_non_empty_scalar(contract$figures[[1L]]$alt_text)
    expect_sha256(contract$figures[[1L]]$source$checksum)
    expect_sha256(contract$figures[[1L]]$provenance$source_hash)
    expect_sha256(contract$provenance$source_hash)

    ids <- vapply(
      contract$themed_report$elements,
      function(element) element$element_id,
      character(1)
    )
    expect_equal(unname(ids), c(
      "table_run_metadata",
      "table_qc_summary",
      "methods_qc_workflow"
    ))
    invisible(lapply(contract$themed_report$elements, function(element) {
      expect_non_empty_scalar(element$caption)
      expect_sha256(element$provenance$source_hash)
    }))

    expect_equal(manifest$schema_version, "flounder.qc_report_manifest.v1")
    expect_equal(manifest$report_id, contract$report_id)
    expect_equal(manifest$element_count, length(ids))
    expect_equal(manifest$figure_count, 1L)
    expect_equal(manifest$contract$sha256, report$manifest$contract$sha256)
    invisible(lapply(manifest$outputs, function(output) {
      expect_true(output$status %in%
                    c("runtime_unavailable", "not_requested",
                      "pending_runtime_binding"))
      expect_true(output$format %in% c("html", "pdf"))
    }))
  })
})
