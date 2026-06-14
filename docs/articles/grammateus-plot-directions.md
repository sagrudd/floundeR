# Grammateus Plot Directions

This vignette shows the two supported plot directions in floundeR’s
Grammateus reporting interface:

- an analyst-owned R `ggplot2` plot is dropped into a Grammateus report
  contract;
- a Grammateus semantic plot specification is rendered by controlled R
  code and then included in a report contract.

The vignette itself is RMarkdown documentation. The report artifacts it
creates are Grammateus-shaped semantic contracts and manifests, and
governed HTML/PDF rendering remains optional through the private
prebuilt Grammateus runtime.

``` r

library(floundeR)
#> floundeR v0.21.19
```

## Prepare Shared QC Evidence

Both directions start with ordinary floundeR QC tables. They are stable
R objects, so they can be inspected interactively, exported to
Synoptikon, or assembled into Grammateus report elements.

``` r

summary_file <- flnDr("sequencing_summary.txt.bz2")
generated_at <- "2026-06-14T15:42:10Z"
run_id <- "example-grammateus-plot-directions"

report_dir <- file.path(tempdir(), "flounder-grammateus-plot-directions")
unlink(report_dir, recursive = TRUE)
dir.create(report_dir, recursive = TRUE)

run_summary <- qc_run_summary(
  summary_file,
  source_id = "example-sequencing-summary"
)
yield_over_time <- qc_yield_over_time(summary_file, resolution_minutes = 15)
run_report_card <- qc_report_card(run_summary)

elements <- grammateus_qc_report_elements(
  qc_summary = run_summary,
  yield_over_time = yield_over_time,
  report_card_findings = run_report_card,
  methods = paste(
    "Run-level QC was generated from the bundled sequencing-summary fixture.",
    "Plot artifacts were produced through the floundeR Grammateus reporting",
    "interface without requiring private runtime assets."
  ),
  provenance = data.frame(
    input_id = "example-sequencing-summary",
    kind = "sequencing_summary",
    uri = summary_file,
    size_bytes = unname(file.info(summary_file)$size),
    stringsAsFactors = FALSE
  ),
  produced_at_utc = generated_at,
  run_id = run_id
)

theme <- grammateus_mnemosyne_theme(
  produced_at_utc = generated_at,
  run_id = run_id
)

run_summary
#> # A tibble: 1 × 20
#>   schema_version        source_id read_count passed_read_count failed_read_count
#>   <chr>                 <chr>          <int>             <int>             <int>
#> 1 flounder.qc_run_summ… example-…      10000              8355              1645
#> # ℹ 15 more variables: pass_fraction <dbl>, total_bases <dbl>,
#> #   passed_bases <dbl>, failed_bases <dbl>, mean_read_length <dbl>,
#> #   median_read_length <dbl>, n50_read_length <dbl>, mean_qscore <dbl>,
#> #   median_qscore <dbl>, channel_count <int>, first_read_start_time <dbl>,
#> #   last_read_start_time <dbl>, run_duration_seconds <dbl>,
#> #   barcode_count <int>, unclassified_read_count <int>
```

## Direction 1: R Plot Dropped Into A Report

Use this path when an analyst has already produced a `ggplot2` plot in R
and needs to hand it to the governed report pipeline with stable
identifiers, captions, alt text, artifact checksums, and source-data
provenance.

``` r

if (requireNamespace("ggplot2", quietly = TRUE)) {
  yield_plot <- ggplot2::ggplot(
    yield_over_time,
    ggplot2::aes(
      x = bin_end_minutes,
      y = cumulative_bases,
      colour = passes_filtering,
      group = passes_filtering
    )
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::labs(
      x = "Elapsed time (minutes)",
      y = "Cumulative bases",
      colour = "Passes filtering"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  r_owned_figure <- grammateus_figure_from_ggplot(
    plot = yield_plot,
    figure_id = "figure_r_owned_yield_over_time",
    caption = "Cumulative sequencing yield over time from an analyst R plot.",
    alt_text = paste(
      "Line plot showing cumulative sequencing bases over elapsed run time,",
      "split by whether reads passed filtering."
    ),
    output_dir = report_dir,
    formats = "png",
    width = 140,
    height = 90,
    units = "mm",
    dpi = 120,
    methods_note = "Generated in R with ggplot2 and wrapped for Grammateus.",
    source_data = yield_over_time,
    produced_at_utc = generated_at,
    run_id = run_id
  )

  r_owned_figure$source[c("format", "path", "checksum")]
} else {
  r_owned_figure <- NULL
  data.frame(
    status = "not_checked",
    reason = "ggplot2 is not installed in this vignette environment."
  )
}
#> $format
#> [1] "png"
#> 
#> $path
#> [1] "/private/var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T/RtmpzSad7k/flounder-grammateus-plot-directions/figure_r_owned_yield_over_time.png"
#> 
#> $checksum
#> [1] "sha256:f9f8b3b82ca147ac024f7c7aaef5de49c1b8a35a4ef1593f9ee5ff8518f5cf00"
```

The figure can be inserted into a report contract without private
runtime assets. `render = "never"` deliberately writes the semantic
contract and manifest only; authorized environments can switch to
`render = "if_available"` or `render = "require"` after validating a
private Grammateus runtime.

``` r

if (!is.null(r_owned_figure)) {
  r_owned_report <- qc_report(
    elements = elements,
    figures = list(r_owned_figure),
    output_dir = file.path(report_dir, "r-owned-report"),
    output = c("html", "pdf"),
    render = "never",
    theme = theme,
    report_id = "report_r_owned_plot_handoff",
    title = "Example QC report with an R-owned plot",
    produced_at_utc = generated_at,
    run_id = run_id
  )

  r_owned_report$manifest[c(
    "schema_version",
    "report_id",
    "element_count",
    "figure_count"
  )]
}
#> $schema_version
#> [1] "flounder.qc_report_manifest.v1"
#> 
#> $report_id
#> [1] "report_r_owned_plot_handoff"
#> 
#> $element_count
#> [1] 5
#> 
#> $figure_count
#> [1] 1
```

## Direction 2: Grammateus Plot Spec Rendered By R

Use this path when a report template, governed workflow, or
Synoptikon-facing process owns the plot definition. floundeR builds a
backend-neutral semantic plot specification from tidy QC data, then the
controlled R backend renders the requested image artifacts in a
deterministic run directory.

``` r

plot_spec <- qc_plot_yield_over_time(
  yield_over_time,
  output = list(
    width_mm = 140,
    height_mm = 90,
    dpi = 120,
    formats = "png"
  ),
  produced_at_utc = generated_at,
  run_id = run_id
)

plot_spec[c("schema_version", "plot_id", "plot_type", "caption")]
#> $schema_version
#> [1] "flounder.grammateus_plot.v1"
#> 
#> $plot_id
#> [1] "plot_yield_over_time"
#> 
#> $plot_type
#> [1] "line"
#> 
#> $caption
#> [1] "Cumulative sequencing yield over time."
plot_spec$mappings
#> $x
#> [1] "bin_end_minutes"
#> 
#> $y
#> [1] "value"
#> 
#> $color
#> [1] "pass_state"
```

The render step needs local R plot dependencies and an `Rscript`
executable. It does not need the private Grammateus HTML/PDF runtime
because this is the controlled R plot backend, not final report
rendering.

``` r

can_render_plot <- requireNamespace("ggplot2", quietly = TRUE) &&
  requireNamespace("jsonlite", quietly = TRUE) &&
  requireNamespace("scales", quietly = TRUE) &&
  requireNamespace("viridisLite", quietly = TRUE) &&
  nzchar(Sys.which("Rscript"))

if (can_render_plot) {
  plot_run <- grammateus_render_plot(
    plot_spec,
    execution = "local_rscript",
    run_root = file.path(report_dir, "plot-runs"),
    required_packages = c("ggplot2", "jsonlite", "scales", "viridisLite")
  )

  spec_rendered_figure <- plot_run$artifacts$png$figure

  list(
    run_dir = plot_run$run_dir,
    script_sha256 = plot_run$script_sha256,
    figure_checksum = spec_rendered_figure$source$checksum
  )
} else {
  plot_run <- NULL
  spec_rendered_figure <- NULL
  data.frame(
    status = "not_checked",
    reason = paste(
      "The controlled R plotting backend needs ggplot2, jsonlite, scales,",
      "viridisLite, and Rscript."
    )
  )
}
#> $run_dir
#> [1] "/private/var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T/RtmpzSad7k/flounder-grammateus-plot-directions/plot-runs/plot_yield_over_time"
#> 
#> $script_sha256
#> [1] "sha256:850d8f12b9dbe02a993eaa1f758fa2071f412f4919d07b9853f9b4627a10191f"
#> 
#> $figure_checksum
#> [1] "sha256:bc4e9910a2e371f661db03519c45c39a7ce1241eaa54d4581dd056166dcfb046"
```

The rendered figure uses the same `flounder_grammateus_figure` wrapper
as an analyst-supplied R plot, so downstream report assembly does not
need a second figure contract.

``` r

if (!is.null(spec_rendered_figure)) {
  spec_owned_report <- qc_report(
    elements = elements,
    figures = list(spec_rendered_figure),
    output_dir = file.path(report_dir, "spec-owned-report"),
    output = "html",
    render = "never",
    theme = theme,
    report_id = "report_spec_owned_plot_render",
    title = "Example QC report with a Grammateus-owned plot spec",
    produced_at_utc = generated_at,
    run_id = run_id
  )

  spec_owned_report$manifest[c(
    "schema_version",
    "report_id",
    "element_count",
    "figure_count"
  )]
}
#> $schema_version
#> [1] "flounder.qc_report_manifest.v1"
#> 
#> $report_id
#> [1] "report_spec_owned_plot_render"
#> 
#> $element_count
#> [1] 5
#> 
#> $figure_count
#> [1] 1
```

## Practical Boundary

The two plot directions are deliberately symmetrical at report assembly
time: both produce semantic figure metadata with stable identifiers,
captions, alt text, checksums, source-data hashes, and run provenance.
The difference is who owns the plot definition:

- R-owned plots are best for interactive analyst work and specialized
  exploratory figures.
- Grammateus-owned plot specs are best for reproducible technical
  reports, governed templates, and Synoptikon report lifecycle handoff.

Neither path downloads data, shells out to a Grammateus CLI, or requires
private Grammateus source code during public package checks.
