# Plan a read-only POD5 collection subdivision

`pod5_subdivide_plan()` asks the in-process Rust extension and curated
`pod5-tools` library API to build a deterministic, read-only plan
describing how a POD5 file, folder, run tree, or manifest could be split
for demonstration, review, or workflow development. This function does
not write POD5 files and does not create output directories.

## Usage

``` r
pod5_subdivide_plan(
  path,
  strategy = c("file-count", "sample-label", "elapsed-time", "read-count"),
  files_per_chunk = 4L,
  seconds_per_chunk = NULL,
  reads_per_chunk = NULL
)
```

## Arguments

- path:

  Character scalar. POD5 file, directory, run tree, or manifest JSON
  file to plan from.

- strategy:

  One of `file-count`, `sample-label`, `elapsed-time`, or `read-count`.

- files_per_chunk:

  Positive integer used by `strategy = "file-count"`.

- seconds_per_chunk:

  Optional positive integer used by `strategy = "elapsed-time"`.

- reads_per_chunk:

  Optional positive integer used by `strategy = "read-count"`.

## Value

A data frame with one row per planned chunk and the columns
`schema_version`, `source`, `strategy`, `target`, `chunk_index`,
`chunk_label`, `file_count`, `total_bytes`, `read_count`,
`relative_paths`, and `warnings`.

## Examples

``` r
signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
run_dir <- file.path(tempdir(), "flounder-subdivide-plan-example")
unlink(run_dir, recursive = TRUE)
dir.create(run_dir, recursive = TRUE)
writeBin(c(signature, as.raw(rep(0, 16)), signature), file.path(run_dir, "a.pod5"))
pod5_subdivide_plan(run_dir, files_per_chunk = 1)
#>   schema_version
#> 1              1
#>                                                                                         source
#> 1 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//RtmpfBxQW1/flounder-subdivide-plan-example
#>     strategy              target chunk_index chunk_label file_count total_bytes
#> 1 file-count 1 file(s) per chunk           1  chunk-0001          1          32
#>   read_count relative_paths warnings
#> 1         NA         a.pod5     <NA>
```
