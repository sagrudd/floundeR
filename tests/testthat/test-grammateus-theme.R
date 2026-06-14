library(floundeR)

test_that("grammateus_mnemosyne_theme defines runtime-free branding metadata", {
  theme <- grammateus_mnemosyne_theme(
    produced_at_utc = "2026-06-14T11:27:09Z",
    run_id = "theme-test-001"
  )

  expect_s3_class(theme, "flounder_grammateus_theme")
  expect_equal(theme$schema_version, "flounder.grammateus_theme.v1")
  expect_equal(theme$brand$name, "Mnemosyne Biosciences")
  expect_equal(theme$template$template_id, "mnemosyne_technical_qc_v1")
  expect_equal(theme$template$runtime_template_ref,
               "grammateus://templates/mnemosyne_technical_qc_v1")
  expect_equal(theme$theme$theme_id, "mnemosyne_biosciences_qc_v1")
  expect_equal(theme$theme$runtime_theme_ref,
               "grammateus://themes/mnemosyne_biosciences_qc_v1")
  expect_equal(theme$theme$profile, "technical_qc")
  expect_equal(theme$theme$palette, "mnemosyne_qc")
  expect_true(theme$style_policy$no_inline_css)
  expect_true(theme$style_policy$no_rmarkdown_css)
  expect_true(theme$style_policy$mnemosyne_branding_required)
  expect_match(theme$theme$status_colours$pass, "^#[0-9A-Fa-f]{6}$")
  expect_match(theme$theme$status_colours$warn, "^#[0-9A-Fa-f]{6}$")
  expect_match(theme$theme$status_colours$fail, "^#[0-9A-Fa-f]{6}$")
  expect_match(theme$theme$status_colours$not_checked, "^#[0-9A-Fa-f]{6}$")
  expect_match(theme$provenance$source_hash, "^sha256:[0-9a-f]{64}$")
  expect_equal(theme$provenance$produced_at_utc, "2026-06-14T11:27:09Z")
  expect_equal(theme$provenance$run_id, "theme-test-001")
})

test_that("grammateus_apply_theme wraps a single semantic report element", {
  element <- grammateus_report_element(
    element_id = "table_qc_summary",
    element_type = "table",
    title = "QC summary",
    caption = "Run-level sequencing summary metrics.",
    data = data.frame(metric = "read_count", value = 3),
    produced_at_utc = "2026-06-14T11:27:09Z"
  )
  theme <- grammateus_mnemosyne_theme(
    produced_at_utc = "2026-06-14T11:27:09Z"
  )

  report <- grammateus_apply_theme(
    element,
    theme = theme,
    produced_at_utc = "2026-06-14T11:27:09Z",
    run_id = "theme-report-001"
  )

  expect_s3_class(report, "flounder_grammateus_themed_report")
  expect_equal(report$schema_version, "flounder.grammateus_themed_report.v1")
  expect_equal(report$element_count, 1L)
  expect_named(report$elements, "table_qc_summary")
  expect_equal(report$theme$brand$name, "Mnemosyne Biosciences")
  expect_s3_class(report$elements$table_qc_summary,
                  "flounder_grammateus_report_element")
  expect_match(report$provenance$source_hash, "^sha256:[0-9a-f]{64}$")
  expect_equal(report$provenance$run_id, "theme-report-001")
})

test_that("grammateus_apply_theme preserves report element bundle names", {
  bundle <- grammateus_qc_report_elements(
    run_metadata = data.frame(run_id = "run-theme"),
    qc_summary = data.frame(metric = "read_count", value = 3),
    methods = "Fixture metrics were summarised by floundeR.",
    produced_at_utc = "2026-06-14T11:27:09Z"
  )

  report <- grammateus_apply_theme(
    bundle,
    produced_at_utc = "2026-06-14T11:27:09Z"
  )

  expect_equal(report$element_count, bundle$element_count)
  expect_named(report$elements, names(bundle$elements))
  expect_equal(report$elements$run_metadata$element_id, "table_run_metadata")
  expect_equal(report$elements$methods$element_type, "methods")
})

test_that("grammateus theme helpers validate profiles and element inputs", {
  expect_error(
    grammateus_mnemosyne_theme(profile = "marketing"),
    "'arg' should be one of"
  )

  expect_error(
    grammateus_mnemosyne_theme(runtime_theme_id = "bad theme"),
    "runtime_theme_id must use"
  )

  expect_error(
    grammateus_apply_theme(
      list(bad = data.frame(x = 1)),
      theme = grammateus_mnemosyne_theme()
    ),
    "all elements must be"
  )

  bad_theme <- grammateus_mnemosyne_theme()
  bad_theme$style_policy$no_inline_css <- FALSE
  expect_error(
    grammateus_apply_theme(
      grammateus_report_element(
        element_id = "table_qc_summary",
        element_type = "table",
        title = "QC summary",
        caption = "Run-level sequencing summary metrics.",
        data = data.frame(metric = "read_count", value = 3)
      ),
      theme = bad_theme
    ),
    "avoid inline CSS"
  )
})
