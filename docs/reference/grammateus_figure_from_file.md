# Build Grammateus figure metadata from image artifacts

`grammateus_figure_from_file()` wraps an existing PNG or SVG file as a
Grammateus-shaped semantic figure. `grammateus_figure_from_ggplot()`
saves a `ggplot2` plot to deterministic PNG and/or SVG artifacts before
wrapping the generated files. These helpers do not require the private
Grammateus runtime; they prepare the open-source R-side handoff object
that later report rendering bindings consume.

## Usage

``` r
grammateus_figure_from_file(
  path,
  figure_id,
  caption,
  alt_text,
  methods_note = NULL,
  layout_hint = c("inline", "full_width", "two_column", "appendix", "panel_group"),
  source_hash = NULL,
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  width_px = NULL,
  height_px = NULL,
  transparent_background = FALSE
)

grammateus_figure_from_ggplot(
  plot,
  figure_id,
  caption,
  alt_text,
  output_dir,
  formats = c("svg", "png"),
  width = 160,
  height = 100,
  units = c("mm", "in", "cm", "px"),
  dpi = 300,
  methods_note = NULL,
  layout_hint = c("inline", "full_width", "two_column", "appendix", "panel_group"),
  source_hash = NULL,
  source_data = NULL,
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  transparent_background = TRUE,
  background = "transparent",
  overwrite = TRUE
)
```

## Arguments

- path:

  Path to an existing PNG or SVG file.

- figure_id:

  Stable lower-snake-case figure identifier. Must start with `figure_`.

- caption:

  Figure caption.

- alt_text:

  Required accessibility text for HTML output.

- methods_note:

  Optional methods or acquisition note.

- layout_hint:

  Grammateus figure layout hint.

- source_hash:

  Optional SHA-256 hash for the source data that produced the figure.
  Plain 64-character hex digests are accepted and normalised to
  `sha256:<hex>`.

- produced_by:

  Producing software or service name.

- producer_version:

  Producing software version.

- produced_at_utc:

  Production timestamp.

- run_id:

  Optional upstream run, workflow, or analysis identifier.

- width_px, height_px:

  Pixel dimensions for an existing file. Required only when SVG
  dimensions cannot be inferred.

- transparent_background:

  Whether transparent background preservation is expected.

- plot:

  A `ggplot2` plot object.

- output_dir:

  Directory where plot artifacts are written.

- formats:

  One or more artifact formats: `png`, `svg`.

- width, height:

  Plot dimensions.

- units:

  Unit for `width` and `height`; one of `in`, `cm`, `mm`, or `px`.

- dpi:

  Output resolution used for PNG output and pixel metadata.

- source_data:

  Optional source data used to compute `source_hash` when wrapping a
  `ggplot2` object. When omitted, the built plot data are used.

- background:

  Background passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

- overwrite:

  Whether existing deterministic artifact paths may be replaced.

## Value

A list with class `flounder_grammateus_figure`. The object follows the
Grammateus `ReportFigure` shape and contains `figure_id`, `caption`,
`alt_text`, `methods_note`, `layout_hint`, `source`, and `provenance`.
When multiple formats are requested from
`grammateus_figure_from_ggplot()`, a `flounder_grammateus_figure_bundle`
containing one figure per format is returned.
