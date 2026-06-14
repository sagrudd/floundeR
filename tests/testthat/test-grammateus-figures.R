test_that("grammateus_figure_from_ggplot writes deterministic PNG metadata", {
  skip_if_not(capabilities("png"))

  output_dir <- tempfile("flounder-grammateus-")
  plot <- ggplot2::ggplot(
    data.frame(minutes = c(0, 1, 2), bases = c(10, 30, 45)),
    ggplot2::aes(minutes, bases)
  ) +
    ggplot2::geom_line()

  figure <- grammateus_figure_from_ggplot(
    plot = plot,
    figure_id = "figure_yield_over_time",
    caption = "Cumulative sequencing yield over time.",
    alt_text = "Line chart showing cumulative base yield over time.",
    output_dir = output_dir,
    formats = "png",
    width = 2,
    height = 1,
    units = "in",
    dpi = 96,
    methods_note = "Generated from sequencing summary records.",
    source_data = data.frame(minutes = c(0, 1, 2), bases = c(10, 30, 45)),
    produced_at_utc = as.POSIXct("2026-06-14 09:00:00", tz = "UTC"),
    run_id = "run-001"
  )

  expect_s3_class(figure, "flounder_grammateus_figure")
  expect_equal(figure$schema_version, "flounder.grammateus_figure.v1")
  expect_equal(figure$figure_id, "figure_yield_over_time")
  expect_equal(figure$source$kind, "image")
  expect_equal(figure$source$format, "png")
  expect_equal(basename(figure$source$path), "figure_yield_over_time.png")
  expect_true(file.exists(figure$source$path))
  expect_match(figure$source$checksum, "^sha256:[0-9a-f]{64}$")
  expect_match(figure$provenance$source_hash, "^sha256:[0-9a-f]{64}$")
  expect_equal(figure$source$width_px, 192L)
  expect_equal(figure$source$height_px, 96L)
  expect_equal(figure$provenance$produced_by, "floundeR")
  expect_equal(figure$provenance$produced_at_utc, "2026-06-14T09:00:00Z")
  expect_equal(figure$provenance$run_id, "run-001")
})

test_that("grammateus_figure_from_file wraps existing SVG metadata", {
  skip_if_not_installed("svglite")

  output_dir <- tempfile("flounder-grammateus-")
  plot <- ggplot2::ggplot(
    data.frame(x = 1:3, y = c(3, 2, 5)),
    ggplot2::aes(x, y)
  ) +
    ggplot2::geom_point()

  bundle <- grammateus_figure_from_ggplot(
    plot = plot,
    figure_id = "figure_quality_distribution",
    caption = "Distribution of read quality.",
    alt_text = "Point plot used as a test figure.",
    output_dir = output_dir,
    formats = c("svg", "png"),
    width = 40,
    height = 30,
    units = "mm",
    dpi = 100,
    source_hash = paste0(rep("a", 64), collapse = ""),
    produced_at_utc = "2026-06-14T09:01:00Z"
  )

  expect_s3_class(bundle, "flounder_grammateus_figure_bundle")
  expect_equal(bundle$formats, c("svg", "png"))
  expect_named(bundle$figures, c("svg", "png"))
  expect_equal(bundle$figures$svg$source$format, "svg")
  expect_equal(bundle$figures$png$source$format, "png")
  expect_match(bundle$figures$svg$source$checksum, "^sha256:[0-9a-f]{64}$")
  expect_equal(
    bundle$figures$svg$provenance$source_hash,
    paste0("sha256:", paste0(rep("a", 64), collapse = ""))
  )

  wrapped <- grammateus_figure_from_file(
    path = bundle$figures$svg$source$path,
    figure_id = "figure_quality_distribution",
    caption = "Distribution of read quality.",
    alt_text = "Point plot used as a test figure.",
    source_hash = bundle$figures$svg$provenance$source_hash,
    produced_at_utc = "2026-06-14T09:02:00Z"
  )

  expect_s3_class(wrapped, "flounder_grammateus_figure")
  expect_equal(wrapped$source$format, "svg")
  expect_equal(wrapped$source$checksum, bundle$figures$svg$source$checksum)
  expect_true(wrapped$source$width_px > 0)
  expect_true(wrapped$source$height_px > 0)
})

test_that("grammateus figure helpers validate Grammateus identifiers and text", {
  plot <- ggplot2::ggplot(
    data.frame(x = 1, y = 1),
    ggplot2::aes(x, y)
  ) +
    ggplot2::geom_point()

  expect_error(
    grammateus_figure_from_ggplot(
      plot = plot,
      figure_id = "plot_bad",
      caption = "Bad identifier.",
      alt_text = "Bad identifier.",
      output_dir = tempfile(),
      formats = "png"
    ),
    "figure_id must start"
  )

  expect_error(
    grammateus_figure_from_ggplot(
      plot = plot,
      figure_id = "figure_bad",
      caption = "",
      alt_text = "Bad caption.",
      output_dir = tempfile(),
      formats = "png"
    ),
    "caption must be"
  )
})
