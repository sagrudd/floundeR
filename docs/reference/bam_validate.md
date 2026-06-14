# Validate BAM structure and selected consistency checks

`bam_validate()` calls Bamana's structural validation library path
in-process. Validation findings are returned as evidence tables even
when Bamana reports validation failure; missing or unreadable inputs
without a validation payload still raise typed Bamana-backed R
conditions.

## Usage

``` r
bam_validate(
  path,
  max_errors = 100L,
  max_warnings = 100L,
  header_only = FALSE,
  records = NULL,
  fail_fast = FALSE,
  include_warnings = TRUE
)
```

## Arguments

- path:

  Character scalar. Local BAM file to validate.

- max_errors:

  Positive integer scalar. Maximum error findings to retain.

- max_warnings:

  Non-negative integer scalar. Maximum warning findings to retain.

- header_only:

  Logical scalar. Validate the header only.

- records:

  Optional positive integer scalar. Validate at most this many records
  after the header.

- fail_fast:

  Logical scalar. Stop at the first error-level finding.

- include_warnings:

  Logical scalar. Include warning-level findings.

## Value

A named list with `status`, `summary`, `findings`, and `error` data
frames. `status$ok` preserves Bamana's command envelope result; a
validation failure with findings therefore returns `ok = FALSE` rather
than raising an R error.

## Examples

``` r
if (FALSE) { # \dontrun{
bam_validate("reads.bam")
} # }
```
