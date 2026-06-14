# Build a library-preparation QC report card

`library_preparation_report_card()` evaluates Porkchop-derived
library-preparation evidence against pass/warn/fail checks for reporting
and review. It consumes R-native evidence tables returned by
[`library_kit_candidates()`](https://sagrudd.github.io/floundeR/reference/library_kit_candidates.md),
[`library_adapter_primer_evidence()`](https://sagrudd.github.io/floundeR/reference/library_kit_candidates.md),
[`library_barcode_evidence()`](https://sagrudd.github.io/floundeR/reference/library_kit_candidates.md),
and
[`library_cdna_primer_evidence()`](https://sagrudd.github.io/floundeR/reference/library_kit_candidates.md);
it does not call Porkchop directly.

## Usage

``` r
library_preparation_report_card(
  kit_candidates = NULL,
  adapter_primer = NULL,
  barcode = NULL,
  cdna = NULL,
  expected_kit_id = NULL,
  thresholds = library_preparation_report_card_thresholds(),
  schema_version = "flounder.library_preparation_report_card.v1"
)
```

## Arguments

- kit_candidates:

  Optional
  [`library_kit_candidates()`](https://sagrudd.github.io/floundeR/reference/library_kit_candidates.md)
  data frame.

- adapter_primer:

  Optional
  [`library_adapter_primer_evidence()`](https://sagrudd.github.io/floundeR/reference/library_kit_candidates.md)
  data frame.

- barcode:

  Optional
  [`library_barcode_evidence()`](https://sagrudd.github.io/floundeR/reference/library_kit_candidates.md)
  data frame.

- cdna:

  Optional
  [`library_cdna_primer_evidence()`](https://sagrudd.github.io/floundeR/reference/library_kit_candidates.md)
  data frame.

- expected_kit_id:

  Optional scalar expected kit id from run metadata or a user-supplied
  review context.

- thresholds:

  Named threshold list. Defaults to
  [`library_preparation_report_card_thresholds()`](https://sagrudd.github.io/floundeR/reference/library_preparation_report_card_thresholds.md).

- schema_version:

  Schema version label to attach to the returned table.

## Value

A tibble using the standard report-card schema: `schema_version`,
`check_id`, `check_label`, `status`, `observed_value`, `warn_threshold`,
`fail_threshold`, `comparator`, and `details`.

## Examples

``` r
kit_candidates <- data.frame(
  kit_id = "SQK-LSK114",
  normalized_score = 0.9,
  support_level = "supported",
  lifecycle_status = "current"
)
library_preparation_report_card(
  kit_candidates = kit_candidates,
  expected_kit_id = "SQK-LSK114"
)
#> # A tibble: 5 × 9
#>   schema_version       check_id check_label status observed_value warn_threshold
#>   <chr>                <chr>    <chr>       <chr>           <dbl>          <dbl>
#> 1 flounder.library_pr… library… Unexpected… pass                0           0   
#> 2 flounder.library_pr… library… Adapter mo… warn               NA           0.2 
#> 3 flounder.library_pr… library… Barcode am… warn               NA           0.05
#> 4 flounder.library_pr… library… cDNA parti… warn               NA           0.2 
#> 5 flounder.library_pr… library… Unsupporte… pass                0           0   
#> # ℹ 3 more variables: fail_threshold <dbl>, comparator <chr>, details <chr>
```
