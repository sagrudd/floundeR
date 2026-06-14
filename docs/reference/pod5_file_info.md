# Inspect local POD5 file metadata

`pod5_file_info()` returns file-level POD5 metadata through the
in-process Rust extension and the curated `pod5-tools` library API.
Until the deep POD5 parser backend is connected, internal fields such as
read count, flow cell ID, sequencing kit, and POD5 version are returned
as missing values, while file size and integrity availability are still
reported.

## Usage

``` r
pod5_file_info(path)
```

## Arguments

- path:

  Character scalar. POD5 file to inspect.

## Value

A one-row data frame with the columns `path`, `size_bytes`,
`flow_cell_id`, `sequencing_kit`, `read_count`, `acquisition_start_utc`,
`duration_seconds`, `pod5_version`, `integrity_status`, and
`integrity_reason`.

## Examples

``` r
signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
path <- file.path(tempdir(), "flounder-file-info-example.pod5")
writeBin(c(signature, as.raw(rep(0, 16)), signature), path)
pod5_file_info(path)
#>                                                                                           path
#> 1 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//RtmpfBxQW1/flounder-file-info-example.pod5
#>   size_bytes flow_cell_id sequencing_kit read_count acquisition_start_utc
#> 1         32         <NA>           <NA>         NA                  <NA>
#>   duration_seconds pod5_version integrity_status
#> 1               NA         <NA>      unavailable
#>                                                             integrity_reason
#> 1 POD5 parser backend not configured; only filesystem metadata was inspected
```
