# Build Grammateus plot specs for nanopore QC report families

These helpers convert floundeR QC tables and curated Bamana/POD5 summary
objects into backend-neutral Grammateus `ReportPlot` specifications.
They do not render plots directly and do not require private Grammateus
runtime assets.

## Usage

``` r
qc_plot_yield_over_time(
  x,
  y = c("cumulative_bases", "cumulative_read_count"),
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  output = NULL
)

qc_plot_quality_distribution(
  x,
  y = c("read_count", "bases"),
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  output = NULL
)

qc_plot_read_length_distribution(
  x,
  y = c("read_count", "bases"),
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  output = NULL
)

qc_plot_flowcell_density(
  x,
  y = c("bases", "read_count", "pass_fraction"),
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  output = NULL
)

qc_plot_barcode_balance(
  x,
  y = c("read_fraction", "bases_fraction", "read_count", "bases"),
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  output = NULL
)

qc_plot_pod5_integrity(
  x,
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  output = NULL
)

qc_plot_bam_mapping_summary(
  x,
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  output = NULL
)

qc_plot_bam_mapq_distribution(
  x,
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  output = NULL
)

qc_plot_bam_flag_summary(
  x,
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  output = NULL
)
```

## Arguments

- x:

  A floundeR QC table, POD5 evidence table, or
  [`bam_summary()`](https://sagrudd.github.io/floundeR/reference/bam_summary.md)
  result appropriate for the plot family.

- y:

  Metric to plot where the plot family supports alternatives.

- produced_by:

  Producing software or service name.

- producer_version:

  Producing software version.

- produced_at_utc:

  Production timestamp.

- run_id:

  Optional upstream run, workflow, or analysis identifier.

- output:

  Grammateus plot output definition.

## Value

A `flounder_grammateus_plot_spec` object.
