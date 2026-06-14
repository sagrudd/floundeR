# Build a versioned POD5 collection manifest

`pod5_manifest()` inventories a local POD5 file, folder, or run tree
through the in-process Rust extension and the curated `pod5-tools`
library API. The returned manifest is intentionally row-oriented for R
and synoptikon handoff: every row carries the manifest schema version
and source path alongside one POD5 entry.

## Usage

``` r
pod5_manifest(path)
```

## Arguments

- path:

  Character scalar. POD5 file or directory to inventory.

## Value

A data frame with one row per POD5 file and the columns
`schema_version`, `source`, `relative_path`, `path`, `size_bytes`,
`verification_status`, and `verification_failed_checks`.

## Examples

``` r
signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
run_dir <- file.path(tempdir(), "flounder-manifest-example")
unlink(run_dir, recursive = TRUE)
dir.create(run_dir, recursive = TRUE)
writeBin(c(signature, as.raw(rep(0, 16)), signature), file.path(run_dir, "a.pod5"))
pod5_manifest(run_dir)
#>   schema_version
#> 1              1
#>                                                                                   source
#> 1 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//RtmpfBxQW1/flounder-manifest-example
#>   relative_path
#> 1        a.pod5
#>                                                                                            path
#> 1 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//RtmpfBxQW1/flounder-manifest-example/a.pod5
#>   size_bytes verification_status verification_failed_checks
#> 1         32          incomplete                          0
```
