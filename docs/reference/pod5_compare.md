# Compare two POD5 collections or manifests

`pod5_compare()` compares two local POD5 files, folders, run trees, or
`pod5-tools` manifest JSON files through the in-process Rust extension.
It is designed for operational handoff checks: missing files and changed
file-size or verification-state differences are returned as rows rather
than rendered text.

## Usage

``` r
pod5_compare(left, right)
```

## Arguments

- left, right:

  Character scalars. POD5 inputs or manifest JSON files to compare.

## Value

A data frame with the columns `status`, `kind`, `relative_path`,
`left_size_bytes`, `right_size_bytes`, `left_verification_status`, and
`right_verification_status`. Matching inputs return one row with
`status = "match"` and `kind = "match"`.

## Examples

``` r
signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
left <- file.path(tempdir(), "flounder-compare-left")
right <- file.path(tempdir(), "flounder-compare-right")
unlink(c(left, right), recursive = TRUE)
dir.create(left, recursive = TRUE)
dir.create(right, recursive = TRUE)
writeBin(c(signature, as.raw(rep(0, 16)), signature), file.path(left, "a.pod5"))
writeBin(c(signature, as.raw(rep(1, 16)), signature), file.path(right, "a.pod5"))
pod5_compare(left, right)
#>   status  kind relative_path left_size_bytes right_size_bytes
#> 1  match match          <NA>              NA               NA
#>   left_verification_status right_verification_status
#> 1                     <NA>                      <NA>
```
