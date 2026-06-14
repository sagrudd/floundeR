# Summarise a local POD5 folder or run tree

`pod5_folder_info()` aggregates file-level POD5 metadata and fast
verification state across a local directory tree through the in-process
Rust extension and the curated `pod5-tools` library API. The current
backend detects duplicate file names and implemented verification
failures while reporting unavailable deep POD5 metadata explicitly in
the `warnings` and `integrity_*` fields.

## Usage

``` r
pod5_folder_info(path)
```

## Arguments

- path:

  Character scalar. Directory containing POD5 files.

## Value

A one-row data frame with the columns `path`, `pod5_file_count`,
`total_bytes`, `total_reads`, `flow_cell_ids`, `sequencing_kits`,
`acquisition_start_utc`, `acquisition_end_utc`, `integrity_status`,
`integrity_reason`, `failed_file_count`, `verification_failed_count`,
`duplicate_file_names`, and `warnings`.

## Examples

``` r
signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
run_dir <- file.path(tempdir(), "flounder-folder-info-example")
unlink(run_dir, recursive = TRUE)
dir.create(run_dir, recursive = TRUE)
writeBin(c(signature, as.raw(rep(0, 16)), signature), file.path(run_dir, "a.pod5"))
pod5_folder_info(run_dir)
#>                                                                                        path
#> 1 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//RtmpfBxQW1/flounder-folder-info-example
#>   pod5_file_count total_bytes total_reads flow_cell_ids sequencing_kits
#> 1               1          32          NA          <NA>            <NA>
#>   acquisition_start_utc acquisition_end_utc integrity_status
#> 1                  <NA>                <NA>      unavailable
#>                                  integrity_reason failed_file_count
#> 1 deep POD5 integrity requires the parser backend                 0
#>   verification_failed_count duplicate_file_names
#> 1                         0                 <NA>
#>                                                                                                                                                                                                                              warnings
#> 1 flow cell metadata unavailable with current POD5 reader backend; sequencing kit metadata unavailable with current POD5 reader backend; acquisition timestamps unavailable for one or more files; temporal gap checks are incomplete
```
