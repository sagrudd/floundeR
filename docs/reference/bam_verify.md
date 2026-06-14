# Verify a BAM file header and container identity

`bam_verify()` calls Bamana's shallow verification library path
in-process. It confirms BGZF container recognition, BAM magic, and
native BAM header parsing, but it does not scan alignment records and
does not prove EOF completeness.

## Usage

``` r
bam_verify(path)
```

## Arguments

- path:

  Character scalar. Local BAM file to verify.

## Value

A one-row data frame with command status, detected format, container,
shallow/deep validation flags, checks performed, confidence, semantic
note, and Bamana error metadata columns.

## Examples

``` r
if (FALSE) { # \dontrun{
bam_verify("reads.bam")
} # }
```
