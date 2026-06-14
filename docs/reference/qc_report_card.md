# Build a run-level QC report card

`qc_report_card()` evaluates a run summary against pass/warn/fail
threshold checks. It is intentionally limited to
sequencing-summary-derived evidence until POD5, BAM/Bamana, and Porkchop
contracts are available.

## Usage

``` r
qc_report_card_thresholds()

qc_report_card(
  x,
  thresholds = qc_report_card_thresholds(),
  barcode_composition = NULL,
  schema_version = "flounder.qc_report_card.v1"
)
```

## Arguments

- x:

  A `SequencingSummary` object, a sequencing-summary file path, a
  normalised sequencing-summary data frame/tibble, or a
  [`qc_run_summary()`](https://sagrudd.github.io/floundeR/reference/qc_run_summary.md)
  one-row tibble.

- thresholds:

  Named threshold list. Defaults to `qc_report_card_thresholds()`.

- barcode_composition:

  Optional
  [`qc_barcode_composition()`](https://sagrudd.github.io/floundeR/reference/qc_tidy_outputs.md)
  tibble. When omitted and `x` is sequencing-summary-like, it is derived
  automatically.

- schema_version:

  Schema version label to attach to the returned table.

## Value

`qc_report_card_thresholds()` returns a named list of warn/fail
thresholds.

`qc_report_card()` returns a tibble with one row per report-card check:

- `schema_version`:

  Report-card schema identifier.

- `check_id`:

  Stable lower-snake-case check identifier.

- `check_label`:

  Human-readable check label.

- `status`:

  One of `pass`, `warn`, or `fail`.

- `observed_value`:

  Numeric value evaluated by the check.

- `warn_threshold`:

  Warning threshold used by the check.

- `fail_threshold`:

  Failure threshold used by the check.

- `comparator`:

  Threshold direction: `minimum` or `maximum`.

- `details`:

  Short explanatory text for report rendering.

## Details

Default thresholds are:

- `pass_fraction_min`:

  warn \< 0.85, fail \< 0.70

- `mean_qscore_min`:

  warn \< 10, fail \< 8

- `n50_read_length_min`:

  warn \< 1000 bases, fail \< 500 bases

- `total_bases_min`:

  warn \< 1,000,000 bases, fail \< 100,000 bases

- `channel_count_min`:

  warn \< 64 channels, fail \< 16 channels

- `unclassified_fraction_max`:

  warn \> 0.10, fail \> 0.25

- `barcode_max_fraction`:

  warn \> 0.60, fail \> 0.80

Missing observed values are reported as `warn` with explanatory detail
text so report consumers can distinguish absent evidence from explicit
failure.

## Examples

``` r
summary_file <- flnDr("sequencing_summary.txt.bz2")
qc_report_card(summary_file)
#> # A tibble: 7 × 9
#>   schema_version       check_id check_label status observed_value warn_threshold
#>   <chr>                <chr>    <chr>       <chr>           <dbl>          <dbl>
#> 1 flounder.qc_report_… pass_fr… Pass read … warn            0.836           0.85
#> 2 flounder.qc_report_… mean_qs… Mean read … warn            9.48           10   
#> 3 flounder.qc_report_… n50_rea… Read lengt… pass        26855            1000   
#> 4 flounder.qc_report_… total_b… Total base… pass    147095512         1000000   
#> 5 flounder.qc_report_… channel… Active cha… pass          502              64   
#> 6 flounder.qc_report_… unclass… Unclassifi… pass            0               0.1 
#> 7 flounder.qc_report_… barcode… Largest ba… warn           NA               0.6 
#> # ℹ 3 more variables: fail_threshold <dbl>, comparator <chr>, details <chr>
```
