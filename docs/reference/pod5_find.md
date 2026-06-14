# Discover local POD5-containing folders

`pod5_find()` recursively searches a local directory tree for folders
that directly contain one or more `.pod5` files. Discovery is
implemented through the in-process Rust extension and the curated
`pod5-tools` library API; it does not invoke an external command-line
process.

## Usage

``` r
pod5_find(path)
```

## Arguments

- path:

  Character scalar. Directory to search.

## Value

A data frame with one row per discovered directory and the columns
`path`, `pod5_file_count`, `total_bytes`, `oldest_modified_utc`, and
`newest_modified_utc`.

## Examples

``` r
run_dir <- file.path(tempdir(), "flounder-pod5-find-example")
unlink(run_dir, recursive = TRUE)
dir.create(file.path(run_dir, "sample-a"), recursive = TRUE)
dir.create(file.path(run_dir, "sample-b"), recursive = TRUE)

writeBin(charToRaw("pod5 placeholder"), file.path(run_dir, "sample-a", "a.pod5"))
writeBin(charToRaw("ignored"), file.path(run_dir, "sample-a", "a.fastq"))
writeBin(charToRaw("pod5 placeholder"), file.path(run_dir, "sample-b", "b.pod5"))

pod5_find(run_dir)
#>                                                                                               path
#> 1 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//RtmpfBxQW1/flounder-pod5-find-example/sample-a
#> 2 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//RtmpfBxQW1/flounder-pod5-find-example/sample-b
#>   pod5_file_count total_bytes                 oldest_modified_utc
#> 1               1          16 2026-06-14T16:56:58.580237684+00:00
#> 2               1          16 2026-06-14T16:56:58.580634605+00:00
#>                   newest_modified_utc
#> 1 2026-06-14T16:56:58.580237684+00:00
#> 2 2026-06-14T16:56:58.580634605+00:00
```
