# Check BAM sort declaration and observed sort evidence

`bam_check_sort()` calls Bamana's in-process Rust sort check and returns
a one-row QC table comparing declared header order with observed record
order.

## Usage

``` r
bam_check_sort(path, sample_records = 1000L, strict = FALSE)
```

## Arguments

- path:

  Character scalar. Local BAM file to inspect.

- sample_records:

  Integer scalar. Use `0` for a full-file scan, or a positive value for
  a bounded evidence scan.

- strict:

  Logical scalar. Ask Bamana to use strict sort interpretation.

## Value

A one-row data frame with declared sort order, observed evidence,
agreement, confidence, and first-violation columns.

## Examples

``` r
if (FALSE) { # \dontrun{
bam_check_sort("reads.bam")
} # }
```
