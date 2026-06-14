# Check canonical BGZF EOF evidence

`bam_check_eof()` reports canonical BGZF EOF-marker evidence through
Bamana's in-process Rust library path. It reports tail completeness only
and does not imply BAM header, record, or aux-field validity.

## Usage

``` r
bam_check_eof(path)
```

## Arguments

- path:

  Character scalar. Local BAM file to inspect.

## Value

A one-row data frame with command status, detected format,
`bgzf_eof_present`, `complete`, semantic note, and Bamana error metadata
columns. Missing EOF is returned as `complete = FALSE` evidence.

## Examples

``` r
if (FALSE) { # \dontrun{
bam_check_eof("reads.bam")
} # }
```
