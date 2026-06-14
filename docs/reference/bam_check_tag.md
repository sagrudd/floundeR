# Check BAM aux-tag evidence for library/QC review

`bam_check_tag()` calls Bamana's in-process Rust aux-tag check. It is
useful for QC report cards that need bounded or full-scan evidence for
tags such as read groups, molecular identifiers, or alignment
annotations.

## Usage

``` r
bam_check_tag(
  path,
  tag,
  sample_records = 1000L,
  full_scan = FALSE,
  require_type = NULL,
  count_hits = FALSE
)
```

## Arguments

- path:

  Character scalar. Local BAM file to inspect.

- tag:

  Character scalar. Two-character BAM aux tag to search for.

- sample_records:

  Integer scalar. Number of records to examine when `full_scan = FALSE`;
  use `0` only with `full_scan = TRUE`.

- full_scan:

  Logical scalar. Scan the full file for absence/presence evidence.

- require_type:

  Optional character scalar. Bamana aux type code to require, for
  example `"Z"` or `"i"`.

- count_hits:

  Logical scalar. Count all matching records during the scan.

## Value

A one-row data frame with tag evidence, mode, confidence, and Bamana
error metadata columns.

## Examples

``` r
if (FALSE) { # \dontrun{
bam_check_tag("reads.bam", "RG", full_scan = TRUE)
} # }
```
