# Build tidy QC data from sequencing-summary inputs

These helpers return stable, schema-versioned tibbles for the first
floundeR QC plot and reporting families. They accept the same input
shapes as
[`qc_run_summary()`](https://sagrudd.github.io/floundeR/reference/qc_run_summary.md):
a `SequencingSummary` object, a sequencing-summary file path, or a
normalised data frame/tibble.

## Usage

``` r
qc_yield_over_time(
  x,
  resolution_minutes = 15,
  schema_version = "flounder.qc_yield_over_time.v1"
)

qc_read_length_distribution(
  x,
  bins = 20,
  schema_version = "flounder.qc_read_length_distribution.v1"
)

qc_quality_distribution(
  x,
  bins = 20,
  schema_version = "flounder.qc_quality_distribution.v1"
)

qc_channel_density(x, schema_version = "flounder.qc_channel_density.v1")

qc_barcode_composition(
  x,
  schema_version = "flounder.qc_barcode_composition.v1"
)
```

## Arguments

- x:

  A `SequencingSummary` object, a sequencing-summary file path, or a
  data frame/tibble with floundeR-normalised sequencing-summary columns.

- resolution_minutes:

  Temporal bin size in minutes.

- schema_version:

  Schema version label to attach to the returned table.

- bins:

  Number of distribution bins to produce.

## Value

`qc_yield_over_time()` returns schema version, bin bounds, pass/fail
state, read count, base yield, and cumulative read/base yield.

`qc_read_length_distribution()` returns schema version, read-length bin
bounds, pass/fail state, read count, and base yield.

`qc_quality_distribution()` returns schema version, Q-score bin bounds,
pass/fail state, read count, and base yield.

`qc_channel_density()` returns schema version, channel, read/base yield,
pass/fail read counts, and pass fraction.

`qc_barcode_composition()` returns schema version, barcode arrangement,
read/base yield, pass/fail read counts, and read/base fractions.

## Examples

``` r
summary_file <- flnDr("sequencing_summary.txt.bz2")
qc_yield_over_time(summary_file)
#> # A tibble: 384 × 10
#>    schema_version            bin_start_seconds bin_end_seconds bin_start_minutes
#>    <chr>                                 <dbl>           <dbl>             <dbl>
#>  1 flounder.qc_yield_over_t…                 0             900                 0
#>  2 flounder.qc_yield_over_t…                 0             900                 0
#>  3 flounder.qc_yield_over_t…               900            1800                15
#>  4 flounder.qc_yield_over_t…               900            1800                15
#>  5 flounder.qc_yield_over_t…              1800            2700                30
#>  6 flounder.qc_yield_over_t…              1800            2700                30
#>  7 flounder.qc_yield_over_t…              2700            3600                45
#>  8 flounder.qc_yield_over_t…              2700            3600                45
#>  9 flounder.qc_yield_over_t…              3600            4500                60
#> 10 flounder.qc_yield_over_t…              3600            4500                60
#> # ℹ 374 more rows
#> # ℹ 6 more variables: bin_end_minutes <dbl>, passes_filtering <lgl>,
#> #   read_count <int>, bases <dbl>, cumulative_read_count <int>,
#> #   cumulative_bases <dbl>
qc_read_length_distribution(summary_file, bins = 10)
#> # A tibble: 16 × 6
#>    schema_version     read_length_bin_start read_length_bin_end passes_filtering
#>    <chr>                              <dbl>               <dbl> <lgl>           
#>  1 flounder.qc_read_…                    0                7799. FALSE           
#>  2 flounder.qc_read_…                    0                7799. TRUE            
#>  3 flounder.qc_read_…                 7799.              15598. FALSE           
#>  4 flounder.qc_read_…                 7799.              15598. TRUE            
#>  5 flounder.qc_read_…                15598.              23396. FALSE           
#>  6 flounder.qc_read_…                15598.              23396. TRUE            
#>  7 flounder.qc_read_…                23396.              31195. FALSE           
#>  8 flounder.qc_read_…                23396.              31195. TRUE            
#>  9 flounder.qc_read_…                31195.              38994  FALSE           
#> 10 flounder.qc_read_…                31195.              38994  TRUE            
#> 11 flounder.qc_read_…                38994               46793. FALSE           
#> 12 flounder.qc_read_…                38994               46793. TRUE            
#> 13 flounder.qc_read_…                46793.              54592. FALSE           
#> 14 flounder.qc_read_…                46793.              54592. TRUE            
#> 15 flounder.qc_read_…                54592.              62390. TRUE            
#> 16 flounder.qc_read_…                70189.              77988  TRUE            
#> # ℹ 2 more variables: read_count <int>, bases <dbl>
qc_quality_distribution(summary_file, bins = 10)
#> # A tibble: 11 × 6
#>    schema_version    qscore_bin_start qscore_bin_end passes_filtering read_count
#>    <chr>                        <dbl>          <dbl> <lgl>                 <int>
#>  1 flounder.qc_qual…             0              1.26 FALSE                   218
#>  2 flounder.qc_qual…             1.26           2.53 FALSE                   401
#>  3 flounder.qc_qual…             2.53           3.79 FALSE                   249
#>  4 flounder.qc_qual…             3.79           5.05 FALSE                   283
#>  5 flounder.qc_qual…             5.05           6.32 FALSE                   336
#>  6 flounder.qc_qual…             6.32           7.58 FALSE                   158
#>  7 flounder.qc_qual…             6.32           7.58 TRUE                    175
#>  8 flounder.qc_qual…             7.58           8.84 TRUE                    708
#>  9 flounder.qc_qual…             8.84          10.1  TRUE                   1270
#> 10 flounder.qc_qual…            10.1           11.4  TRUE                   3555
#> 11 flounder.qc_qual…            11.4           12.6  TRUE                   2647
#> # ℹ 1 more variable: bases <dbl>
qc_channel_density(summary_file)
#> # A tibble: 502 × 7
#>    schema_version  channel read_count  bases passed_read_count failed_read_count
#>    <chr>             <int>      <int>  <dbl>             <int>             <int>
#>  1 flounder.qc_ch…       1         15 253786                12                 3
#>  2 flounder.qc_ch…       2         20 279012                19                 1
#>  3 flounder.qc_ch…       3         17 136937                14                 3
#>  4 flounder.qc_ch…       4         12 244863                11                 1
#>  5 flounder.qc_ch…       5         12 202759                11                 1
#>  6 flounder.qc_ch…       6         13 272392                12                 1
#>  7 flounder.qc_ch…       7         17 208594                12                 5
#>  8 flounder.qc_ch…       8         26 418172                24                 2
#>  9 flounder.qc_ch…       9         25 334420                22                 3
#> 10 flounder.qc_ch…      10         31 399163                29                 2
#> # ℹ 492 more rows
#> # ℹ 1 more variable: pass_fraction <dbl>
qc_barcode_composition(summary_file)
#> # A tibble: 1 × 8
#>   schema_version         barcode_arrangement read_count  bases passed_read_count
#>   <chr>                  <chr>                    <int>  <dbl>             <int>
#> 1 flounder.qc_barcode_c… unclassified             10000 1.47e8              8355
#> # ℹ 3 more variables: failed_read_count <int>, read_fraction <dbl>,
#> #   bases_fraction <dbl>
```
