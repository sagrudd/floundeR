test_that("grammateus_plot_spec builds a ReportPlot-shaped line spec", {
  spec <- grammateus_plot_spec(
    plot_id = "plot_yield_over_time",
    plot_type = "line",
    data = data.frame(minutes = c(0, 1, 2), bases = c(10, 30, 45)),
    mappings = list(x = "minutes", y = "bases"),
    axes = list(
      x = list(label = "Elapsed time", unit = "minutes"),
      y = list(label = "Yield", unit = "bases")
    ),
    caption = "Cumulative sequencing yield over time.",
    output = list(width_mm = 80, height_mm = 50, dpi = 100, formats = "png"),
    produced_at_utc = as.POSIXct("2026-06-14 10:00:00", tz = "UTC")
  )

  expect_s3_class(spec, "flounder_grammateus_plot_spec")
  expect_equal(spec$schema_version, "flounder.grammateus_plot.v1")
  expect_equal(spec$plot_id, "plot_yield_over_time")
  expect_equal(spec$plot_type, "line")
  expect_equal(spec$data$kind, "inline_tidy")
  expect_length(spec$data$records, 3)
  expect_equal(spec$mappings$x, "minutes")
  expect_equal(spec$axes$x$scale, "linear")
  expect_equal(spec$output$formats, "png")
  expect_match(spec$provenance$source_hash, "^sha256:[0-9a-f]{64}$")
  expect_equal(spec$provenance$produced_at_utc, "2026-06-14T10:00:00Z")
})

test_that("grammateus_render_plot executes the controlled local R backend", {
  skip_if_not(capabilities("png"))
  skip_if_not_installed("svglite")
  skip_if_not_installed("viridisLite")

  run_root <- tempfile("flounder-grammateus-runs-")
  spec <- grammateus_plot_spec(
    plot_id = "plot_quality_distribution",
    plot_type = "bar",
    data = data.frame(qscore_bin = c("10-15", "15-20"), reads = c(4, 9)),
    mappings = list(x = "qscore_bin", y = "reads"),
    axes = list(
      x = list(label = "Mean Q-score bin"),
      y = list(label = "Reads")
    ),
    caption = "Distribution of mean read quality.",
    output = list(width_mm = 50, height_mm = 35, dpi = 100,
                  formats = c("png", "svg")),
    bar_value_semantics = "counts",
    produced_at_utc = "2026-06-14T10:01:00Z",
    run_id = "run-plot-001"
  )

  run <- grammateus_render_plot(
    spec,
    execution = "local_rscript",
    run_root = run_root
  )

  expect_s3_class(run, "flounder_grammateus_plot_run")
  expect_equal(run$schema_version, "flounder.grammateus_plot_run.v1")
  expect_equal(basename(run$run_dir), "plot_quality_distribution")
  expect_true(file.exists(run$plot_spec_path))
  expect_true(file.exists(run$plot_data_path))
  expect_true(file.exists(run$script_path))
  expect_true(file.exists(run$backend_session_path))
  expect_match(run$plot_spec_sha256, "^sha256:[0-9a-f]{64}$")
  expect_match(run$plot_data_sha256, "^sha256:[0-9a-f]{64}$")
  expect_match(run$script_sha256, "^sha256:[0-9a-f]{64}$")
  expect_named(run$artifacts, c("png", "svg"))
  expect_true(file.exists(run$artifacts$png$path))
  expect_true(file.exists(run$artifacts$svg$path))
  expect_match(run$artifacts$png$sha256, "^sha256:[0-9a-f]{64}$")
  expect_gt(run$artifacts$png$byte_len, 0)
  expect_s3_class(run$artifacts$png$figure, "flounder_grammateus_figure")
  expect_equal(run$artifacts$png$figure$figure_id, "figure_quality_distribution")
  expect_equal(run$artifacts$png$figure$provenance$run_id, "run-plot-001")
})

test_that("grammateus plot specs validate semantic contract requirements", {
  expect_error(
    grammateus_plot_spec(
      plot_id = "figure_bad",
      plot_type = "line",
      data = data.frame(x = 1, y = 1),
      mappings = list(x = "x", y = "y"),
      axes = list(x = list(label = "X"), y = list(label = "Y")),
      caption = "Bad identifier."
    ),
    "plot_id must start"
  )

  expect_error(
    grammateus_plot_spec(
      plot_id = "plot_bad_bar",
      plot_type = "bar",
      data = data.frame(x = "a", y = 1),
      mappings = list(x = "x", y = "y"),
      axes = list(x = list(label = "X"), y = list(label = "Y")),
      caption = "Bad bar."
    ),
    "bar and stacked_bar plots must declare"
  )

  expect_error(
    grammateus_plot_spec(
      plot_id = "plot_bad_numeric",
      plot_type = "line",
      data = data.frame(x = "a", y = 1),
      mappings = list(x = "x", y = "y"),
      axes = list(x = list(label = "X"), y = list(label = "Y")),
      caption = "Bad numeric field."
    ),
    "must be numeric"
  )
})

test_that("grammateus_render_plot reports backend failures as typed conditions", {
  spec <- grammateus_plot_spec(
    plot_id = "plot_backend_failure",
    plot_type = "scatter",
    data = data.frame(x = 1, y = 1),
    mappings = list(x = "x", y = "y"),
    axes = list(x = list(label = "X"), y = list(label = "Y")),
    caption = "Backend failure test.",
    output = list(width_mm = 20, height_mm = 20, dpi = 100, formats = "png")
  )

  expect_error(
    grammateus_render_plot(
      spec,
      run_root = tempfile("flounder-grammateus-fail-"),
      rscript_path = "definitely-not-an-rscript"
    ),
    class = "flounder_grammateus_backend_error"
  )
})
