# Check whether BAM mapping evidence is present

`bam_check_map()` calls Bamana's in-process Rust mapping check. It
reports whether mapped reads are evidenced by a supported index or by a
bounded/full record scan, with reference-level evidence for QC review.

## Usage

``` r
bam_check_map(path, sample_records = 1000L, prefer_index = TRUE)
```

## Arguments

- path:

  Character scalar. Local BAM file to inspect.

- sample_records:

  Integer scalar. Use `0` for a full-file scan, or a positive value for
  a bounded evidence scan.

- prefer_index:

  Logical scalar. Prefer supported index-derived evidence where
  possible.

## Value

A named list with `status`, `index`, `summary`, and `references` data
frames.

## Examples

``` r
if (FALSE) { # \dontrun{
bam_check_map("reads.bam", sample_records = 1000)
} # }
```
