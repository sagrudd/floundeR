# Summarise a BAM file for nanopore QC

`bam_summary()` calls the curated `../bamana` Rust library in-process
and returns R-native tables for QC, review, provenance, and report
assembly. The binding is intentionally focused on read-only operational
summary evidence: it does not attempt to expose all Bamana commands.

## Usage

``` r
bam_summary(
  path,
  sample_records = 0L,
  prefer_index = FALSE,
  include_mapq_hist = FALSE,
  include_flags = TRUE,
  allow_incomplete = FALSE
)
```

## Arguments

- path:

  Character scalar. Local BAM file to summarise.

- sample_records:

  Integer scalar. Use `0` for a full-file scan, or a positive value for
  a bounded evidence scan.

- prefer_index:

  Logical scalar. Prefer index-derived counts where Bamana can use a
  supported sidecar index.

- include_mapq_hist:

  Logical scalar. Include a MAPQ histogram table.

- include_flags:

  Logical scalar. Include flag-category evidence.

- allow_incomplete:

  Logical scalar. Allow Bamana to return partial evidence for incomplete
  BGZF/BAM input when supported by the backend.

## Value

A named list with data frames including `status`, `evidence`, `header`,
`counts`, `fractions`, `fractions_observed`, `mapq`, `mapping`,
`anomalies`, `flag_categories`, `references`, `index_derived`, and
`mapq_histogram`.

## Examples

``` r
if (FALSE) { # \dontrun{
bam_summary("reads.bam")
} # }
```
