# Summarise a nanopore run for QC handoff

`qc_run_summary()` converts a normalised sequencing-summary table into a
stable one-row QC contract suitable for report cards, Grammateus
reports, and future Synoptikon handoff.

## Usage

``` r
qc_run_summary(
  x,
  source_id = NULL,
  schema_version = "flounder.qc_run_summary.v1"
)
```

## Arguments

- x:

  A `SequencingSummary` object, a sequencing-summary file path, or a
  data frame/tibble with floundeR-normalised sequencing-summary columns.

- source_id:

  Optional run, sample, file, or source identifier. When `x` is a path
  and `source_id` is not supplied, the path is used.

- schema_version:

  Schema version label to attach to the returned table.

## Value

A one-row tibble with schema version, source identity, read/yield
totals, pass/fail counts, length and Q-score summaries, channel count,
observed time bounds, and barcode counts.

## Examples

``` r
summary_file <- flnDr("sequencing_summary.txt.bz2")
qc_run_summary(summary_file)
#> # A tibble: 1 × 20
#>   schema_version        source_id read_count passed_read_count failed_read_count
#>   <chr>                 <chr>          <int>             <int>             <int>
#> 1 flounder.qc_run_summ… /private…      10000              8355              1645
#> # ℹ 15 more variables: pass_fraction <dbl>, total_bases <dbl>,
#> #   passed_bases <dbl>, failed_bases <dbl>, mean_read_length <dbl>,
#> #   median_read_length <dbl>, n50_read_length <dbl>, mean_qscore <dbl>,
#> #   median_qscore <dbl>, channel_count <int>, first_read_start_time <dbl>,
#> #   last_read_start_time <dbl>, run_duration_seconds <dbl>,
#> #   barcode_count <int>, unclassified_read_count <int>
```
