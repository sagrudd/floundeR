# Check BAM index evidence for QC report cards

`bam_check_index()` calls Bamana's in-process Rust index inspection path
and returns stable R tables for index presence, support level,
freshness, and candidate sidecars. It is read-only and does not create
or repair indexes.

## Usage

``` r
bam_check_index(path, require = FALSE, prefer_csi = FALSE)
```

## Arguments

- path:

  Character scalar. Local BAM file to inspect.

- require:

  Logical scalar. Treat a missing index as a Bamana command failure
  while still returning structured index evidence when available.

- prefer_csi:

  Logical scalar. Prefer CSI over BAI when both sidecars are present and
  usable.

## Value

A named list with `status`, `index`, `candidates`, and `error` data
frames.

## Examples

``` r
if (FALSE) { # \dontrun{
bam_check_index("reads.bam")
} # }
```
