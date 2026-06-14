# Build a BAM QC report card

`bam_qc_report_card()` evaluates Bamana-derived BAM evidence against
pass/warn/fail checks for alignment QC and reporting. It consumes the
R-native objects returned by
[`bam_summary()`](https://sagrudd.github.io/floundeR/reference/bam_summary.md),
[`bam_check_index()`](https://sagrudd.github.io/floundeR/reference/bam_check_index.md),
[`bam_check_sort()`](https://sagrudd.github.io/floundeR/reference/bam_check_sort.md),
[`bam_check_tag()`](https://sagrudd.github.io/floundeR/reference/bam_check_tag.md),
[`bam_check_eof()`](https://sagrudd.github.io/floundeR/reference/bam_check_eof.md),
and
[`bam_validate()`](https://sagrudd.github.io/floundeR/reference/bam_validate.md);
it does not call Bamana directly.

## Usage

``` r
bam_qc_report_card(
  summary = NULL,
  index = NULL,
  mapping = NULL,
  sorting = NULL,
  tags = NULL,
  eof = NULL,
  validation = NULL,
  provenance_anomalies = NULL,
  expected_tags = character(),
  thresholds = bam_qc_report_card_thresholds(),
  schema_version = "flounder.bam_qc_report_card.v1"
)
```

## Arguments

- summary:

  Optional
  [`bam_summary()`](https://sagrudd.github.io/floundeR/reference/bam_summary.md)
  result.

- index:

  Optional
  [`bam_check_index()`](https://sagrudd.github.io/floundeR/reference/bam_check_index.md)
  result.

- mapping:

  Optional
  [`bam_check_map()`](https://sagrudd.github.io/floundeR/reference/bam_check_map.md)
  result used as a fallback for mapping-fraction evidence when `summary`
  is absent.

- sorting:

  Optional
  [`bam_check_sort()`](https://sagrudd.github.io/floundeR/reference/bam_check_sort.md)
  result.

- tags:

  Optional
  [`bam_check_tag()`](https://sagrudd.github.io/floundeR/reference/bam_check_tag.md)
  data frame or list of such data frames.

- eof:

  Optional
  [`bam_check_eof()`](https://sagrudd.github.io/floundeR/reference/bam_check_eof.md)
  data frame.

- validation:

  Optional
  [`bam_validate()`](https://sagrudd.github.io/floundeR/reference/bam_validate.md)
  result.

- provenance_anomalies:

  Optional numeric anomaly count or data frame of provenance findings
  from a future Bamana provenance/forensic surface.

- expected_tags:

  Character vector of BAM aux tags expected to be present.

- thresholds:

  Named threshold list. Defaults to
  [`bam_qc_report_card_thresholds()`](https://sagrudd.github.io/floundeR/reference/bam_qc_report_card_thresholds.md).

- schema_version:

  Schema version label to attach to the returned table.

## Value

A tibble using the standard report-card schema: `schema_version`,
`check_id`, `check_label`, `status`, `observed_value`, `warn_threshold`,
`fail_threshold`, `comparator`, and `details`.

## Examples

``` r
if (FALSE) { # \dontrun{
bam <- bam_summary("reads.bam")
bam_qc_report_card(summary = bam)
} # }
```
