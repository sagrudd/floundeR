# Grammateus Technical Report

This vignette shows how floundeR prepares a Grammateus-backed technical
QC report. The package vignette is rendered with RMarkdown, but the
report artifact it creates is not an RMarkdown report: it is a
Grammateus semantic report contract plus a manifest, provenance, figure
metadata, and Mnemosyne Biosciences theme references.

The public open-source package can build these contracts without private
Grammateus runtime assets. Governed HTML/PDF rendering is enabled only
when an authorized prebuilt Grammateus runtime is installed and
validated.

``` r

library(floundeR)
#> floundeR v0.21.19
```

## Prepare QC Evidence

Start with run-level QC tables that are stable R objects and can also
feed Synoptikon handoff.

``` r

summary_file <- flnDr("sequencing_summary.txt.bz2")
generated_at <- "2026-06-14T15:27:10Z"
run_id <- "example-grammateus-report"

run_summary <- qc_run_summary(
  summary_file,
  source_id = "example-sequencing-summary"
)
yield_over_time <- qc_yield_over_time(summary_file, resolution_minutes = 15)
quality_distribution <- qc_quality_distribution(summary_file, bins = 5)
flowcell_density <- qc_channel_density(summary_file)
barcode_balance <- qc_barcode_composition(summary_file)
run_report_card <- qc_report_card(
  run_summary,
  barcode_composition = barcode_balance
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
run_report_card[, c("check_id", "status", "observed_value", "details")]
#> # A tibble: 7 × 4
#>   check_id              status observed_value details                           
#>   <chr>                 <chr>           <dbl> <chr>                             
#> 1 pass_fraction         warn            0.836 Observed value crossed the warnin…
#> 2 mean_qscore           warn            9.48  Observed value crossed the warnin…
#> 3 n50_read_length       pass        26855     Observed value is within the pass…
#> 4 total_bases           pass    147095512     Observed value is within the pass…
#> 5 channel_count         pass          502     Observed value is within the pass…
#> 6 unclassified_fraction pass            0     Observed value is within the pass…
#> 7 barcode_max_fraction  warn           NA     Observed value is missing; review…
```

POD5, BAM, and library-preparation evidence remain sectioned report
inputs. This small vignette records POD5 open-data provenance and marks
missing BAM/library evidence explicitly rather than inventing bundled
files.

``` r

pod5_example <- ont_zymo_pod5_example_objects(role = "pass")

pod5_integrity <- data.frame(
  schema_version = "flounder.pod5_open_data_example.v1",
  status = "not_downloaded",
  role = pod5_example$role,
  s3_uri = pod5_example$s3_uri,
  size_bytes = pod5_example$size,
  intended_use = pod5_example$intended_use,
  stringsAsFactors = FALSE
)

bam_alignment_summary <- data.frame(
  schema_version = "flounder.bam_alignment_summary.v1",
  status = "not_checked",
  detail = "No BAM file is bundled with this vignette fixture.",
  stringsAsFactors = FALSE
)

library_preparation <- data.frame(
  schema_version = "flounder.library_preparation_qc.v1",
  status = "not_checked",
  detail = paste(
    "No read-level library-preparation screen is bundled with this vignette.",
    "Porkchop-derived adapter, primer, barcode, kit, and cDNA evidence should",
    "populate this table when available."
  ),
  stringsAsFactors = FALSE
)
```

## Create A Governed Figure Wrapper

Figures handed to Grammateus need stable identifiers, captions, alt
text, artifact checksums, and provenance. This example wraps a tiny SVG
fixture so the vignette remains dependency-light; the R-to-Grammateus
plot handoff vignette covers `ggplot2` figure generation in more detail.

``` r

report_dir <- file.path(tempdir(), "flounder-grammateus-report")
unlink(report_dir, recursive = TRUE)
dir.create(report_dir, recursive = TRUE)

yield_svg <- file.path(report_dir, "yield-over-time.svg")
writeLines(
  c(
    "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"240\" height=\"120\">",
    "  <rect width=\"240\" height=\"120\" fill=\"#ffffff\"/>",
    "  <polyline points=\"20,95 75,72 130,48 185,34 220,24\"",
    "    fill=\"none\" stroke=\"#237A57\" stroke-width=\"5\"/>",
    "  <circle cx=\"220\" cy=\"24\" r=\"6\" fill=\"#237A57\"/>",
    "</svg>"
  ),
  yield_svg,
  useBytes = TRUE
)

yield_figure <- grammateus_figure_from_file(
  path = yield_svg,
  figure_id = "figure_yield_over_time",
  caption = "Cumulative sequencing yield over time.",
  alt_text = paste(
    "Line chart icon showing cumulative sequencing yield increasing over",
    "elapsed run time."
  ),
  methods_note = "Wrapped as a deterministic SVG artifact for vignette use.",
  width_px = 240,
  height_px = 120,
  produced_at_utc = generated_at,
  run_id = run_id
)

yield_figure$source[c("format", "checksum", "width_px", "height_px")]
#> $format
#> [1] "svg"
#> 
#> $checksum
#> [1] "sha256:ef04d165342ba4276f0b256f584b039ded405917087fba85b152c49c7e0da015"
#> 
#> $width_px
#> [1] 240
#> 
#> $height_px
#> [1] 120
```

## Plot Report-Card Findings

Report-card findings are ordinary tidy rows, so they can be visualised
before they are handed into Grammateus. Two useful examples are a
compact status-count summary for dashboards and a check-level status
matrix for technical reports. Both remain governed figures once wrapped
with
[`grammateus_figure_from_ggplot()`](https://sagrudd.github.io/floundeR/reference/grammateus_figure_from_file.md).

``` r

if (requireNamespace("ggplot2", quietly = TRUE)) {
  status_levels <- c("pass", "warn", "fail", "not_checked")
  status_colours <- c(
    pass = "#237A57",
    warn = "#C88719",
    fail = "#B33A3A",
    not_checked = "#68717A"
  )

  status_counts <- as.data.frame(
    table(
      factor(run_report_card$status, levels = status_levels),
      useNA = "no"
    ),
    stringsAsFactors = FALSE
  )
  names(status_counts) <- c("status", "check_count")
  status_counts <- status_counts[status_counts$check_count > 0, , drop = FALSE]

  status_count_plot <- ggplot2::ggplot(
    status_counts,
    ggplot2::aes(x = status, y = check_count, fill = status)
  ) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::scale_fill_manual(values = status_colours, drop = FALSE) +
    ggplot2::labs(
      x = "Report-card status",
      y = "Configured checks"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "none")

  status_count_figure <- grammateus_figure_from_ggplot(
    plot = status_count_plot,
    figure_id = "figure_report_card_status_counts",
    caption = "Report-card checks by pass, warning, and failure status.",
    alt_text = paste(
      "Bar chart counting report-card checks by pass, warning, failure,",
      "or not-checked status."
    ),
    output_dir = report_dir,
    formats = "png",
    width = 120,
    height = 75,
    units = "mm",
    dpi = 120,
    methods_note = "Generated from qc_report_card() findings.",
    source_data = run_report_card,
    produced_at_utc = generated_at,
    run_id = run_id
  )

  check_status_data <- run_report_card
  check_status_data$check_label <- factor(
    check_status_data$check_label,
    levels = rev(check_status_data$check_label)
  )

  check_status_plot <- ggplot2::ggplot(
    check_status_data,
    ggplot2::aes(x = "QC checks", y = check_label, fill = status)
  ) +
    ggplot2::geom_tile(width = 0.65, height = 0.75) +
    ggplot2::scale_fill_manual(values = status_colours, drop = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = "Status"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  check_status_figure <- grammateus_figure_from_ggplot(
    plot = check_status_plot,
    figure_id = "figure_report_card_check_status",
    caption = "Per-check report-card status for configured run-level QC checks.",
    alt_text = paste(
      "Tile chart listing each report-card check with colour-coded pass,",
      "warning, failure, or not-checked status."
    ),
    output_dir = report_dir,
    formats = "png",
    width = 130,
    height = 90,
    units = "mm",
    dpi = 120,
    methods_note = "Generated from qc_report_card() findings.",
    source_data = run_report_card,
    produced_at_utc = generated_at,
    run_id = run_id
  )

  report_card_figures <- list(status_count_figure, check_status_figure)
  vapply(report_card_figures, function(x) x$source$checksum, character(1))
} else {
  report_card_figures <- list()
  data.frame(
    status = "not_checked",
    reason = "ggplot2 is not installed in this vignette environment."
  )
}
#> [1] "sha256:784bb1d9d8989e2d8ab3b52517fd93c20e8f9e187a5c8706ad24c091c0c0126c"
#> [2] "sha256:b183fddfb0d8b81a615c3d5173618ea05994dced10bf0dcd41a31079a5d00fb2"
```

## Assemble Semantic Report Elements

[`grammateus_qc_report_elements()`](https://sagrudd.github.io/floundeR/reference/grammateus_report_element.md)
maps floundeR QC tables into named semantic technical-report elements.
The Mnemosyne theme descriptor records the private runtime
template/theme references without bundling private assets.

``` r

run_metadata <- data.frame(
  run_id = run_id,
  sample_id = "example-sample",
  project_id = "flounder-revival",
  pod5_source = pod5_example$s3_uri,
  stringsAsFactors = FALSE
)

limitations <- data.frame(
  limitation_id = "example_pod5_not_downloaded",
  severity = "info",
  section = "pod5",
  message = "The selected ONT POD5 object is recorded as provenance only.",
  stringsAsFactors = FALSE
)

input_provenance <- data.frame(
  input_id = "example-sequencing-summary",
  kind = "sequencing_summary",
  uri = summary_file,
  size_bytes = unname(file.info(summary_file)$size),
  stringsAsFactors = FALSE
)

elements <- grammateus_qc_report_elements(
  run_metadata = run_metadata,
  qc_summary = run_summary,
  flowcell_density = flowcell_density,
  yield_over_time = yield_over_time,
  quality_distribution = quality_distribution,
  barcode_balance = barcode_balance,
  pod5_integrity = pod5_integrity,
  bam_alignment_summary = bam_alignment_summary,
  library_preparation = library_preparation,
  report_card_findings = run_report_card,
  methods = paste(
    "Run-level metrics were generated from a floundeR sequencing-summary",
    "fixture. POD5 source metadata references the canonical ONT Zymo fecal",
    "pass example without downloading raw POD5 bytes."
  ),
  limitations = limitations,
  provenance = input_provenance,
  produced_at_utc = generated_at,
  run_id = run_id
)

theme <- grammateus_mnemosyne_theme(
  produced_at_utc = generated_at,
  run_id = run_id
)

elements$element_count
#> [1] 13
theme$brand$name
#> [1] "Mnemosyne Biosciences"
theme$template$runtime_template_ref
#> [1] "grammateus://templates/mnemosyne_technical_qc_v1"
```

## Write The Report Contract

[`qc_report()`](https://sagrudd.github.io/floundeR/reference/qc_report.md)
writes a Grammateus-shaped contract and manifest. The example uses
`render = "never"` so public vignette builds never require private
runtime assets. Use `render = "if_available"` or `render = "require"` in
authorized runtime environments.

``` r

report <- qc_report(
  elements = elements,
  figures = c(list(yield_figure), report_card_figures),
  output_dir = report_dir,
  output = c("html", "pdf"),
  render = "never",
  theme = theme,
  report_id = "report_example_nanopore_qc",
  title = "Example nanopore sequencing QC report",
  produced_at_utc = generated_at,
  run_id = run_id
)

basename(report$contract_path)
#> [1] "report_example_nanopore_qc-contract.json"
basename(report$manifest_path)
#> [1] "report_example_nanopore_qc-manifest.json"
report$manifest[c("schema_version", "report_id", "element_count", "figure_count")]
#> $schema_version
#> [1] "flounder.qc_report_manifest.v1"
#> 
#> $report_id
#> [1] "report_example_nanopore_qc"
#> 
#> $element_count
#> [1] 13
#> 
#> $figure_count
#> [1] 3
vapply(report$outputs, `[[`, character(1), "status")
#> [1] "not_requested" "not_requested"
```

The manifest records the contract checksum, output status, theme
references, and provenance needed for governed report lifecycle
handling.

``` r

manifest <- jsonlite::fromJSON(report$manifest_path, simplifyVector = FALSE)

manifest$contract[c("sha256", "byte_len")]
#> $sha256
#> [1] "sha256:0af39102eb51e62701f3ed4192a61082ada37ec54070939fb921dc05d80a3d8a"
#> 
#> $byte_len
#> [1] 353757
manifest$theme
#> $brand
#> [1] "Mnemosyne Biosciences"
#> 
#> $template_id
#> [1] "mnemosyne_technical_qc_v1"
#> 
#> $runtime_template_ref
#> [1] "grammateus://templates/mnemosyne_technical_qc_v1"
#> 
#> $theme_id
#> [1] "mnemosyne_biosciences_qc_v1"
#> 
#> $runtime_theme_ref
#> [1] "grammateus://themes/mnemosyne_biosciences_qc_v1"
manifest$runtime$status
#> NULL
```

## Runtime-Backed Rendering

Authorized environments can request runtime-backed HTML/PDF rendering
after the prebuilt Grammateus runtime has been installed and validated.

``` r

validation <- grammateus_runtime_validate()

if (isTRUE(validation$ok)) {
  rendered <- qc_report(
    elements = elements,
    figures = c(list(yield_figure), report_card_figures),
    output_dir = report_dir,
    output = c("html", "pdf"),
    render = "if_available",
    report_id = "report_example_nanopore_qc",
    title = "Example nanopore sequencing QC report",
    produced_at_utc = generated_at,
    run_id = run_id
  )
}
```

Use `render = "require"` only when a missing or invalid private runtime
should be treated as a hard failure.
