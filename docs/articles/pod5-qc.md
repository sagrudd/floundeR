# POD5 QC

This vignette shows the package-level POD5 QC workflow in floundeR. It
is a package tutorial, not a production report renderer. Governed
HTML/PDF report artifacts remain a Grammateus responsibility.

The default workflow is network-free and does not download or commit
POD5 bytes. Rust-backed POD5 inspection is demonstrated only when the
compiled extension is available in the local R installation.

``` r

library(floundeR)
#> floundeR v0.21.19
rust_available <- flounder_rust_available()
```

## Canonical Open-Data Source

floundeR uses a fixed ONT open-data POD5 source for opt-in real-data
examples. This is the canonical real-data source for package
documentation, integration checks, report demonstrations, and Synoptikon
examples: `s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/`.

The routine pass example is `PAU85136_pass_279c9095_68316534_8289.pod5`.
The fail object is reserved for fail-state examples and should not be
used as the default tutorial input.

``` r

pod5_dataset <- ont_zymo_pod5_dataset()
pod5_examples <- ont_zymo_pod5_example_objects()

pod5_dataset
#>                           dataset_name        bucket    region
#> 1 ont_zymo_fecal_2025_05_pau85136_pod5 ont-open-data eu-west-1
#>                                  prefix
#> 1 zymo_fecal_2025.05/raw/PAU85136/pod5/
#>                                                     s3_uri total_pod5_objects
#> 1 s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/               9107
#>   pass_pod5_objects fail_pod5_objects  total_bytes         verified_utc
#> 1              8290               817 2.676974e+12 2026-06-13T00:00:00Z
pod5_examples
#>   role state                                 file_name        bucket    region
#> 1 pass  pass PAU85136_pass_279c9095_68316534_8289.pod5 ont-open-data eu-west-1
#> 2 fail  fail    PAU85136_fail_279c9095_68316534_0.pod5 ont-open-data eu-west-1
#>                                                                              key
#> 1 zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_pass_279c9095_68316534_8289.pod5
#> 2    zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_fail_279c9095_68316534_0.pod5
#>                                                                                              s3_uri
#> 1 s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_pass_279c9095_68316534_8289.pod5
#> 2    s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_fail_279c9095_68316534_0.pod5
#>        size        last_modified_utc
#> 1  47077200 2025-05-19T23:24:03.000Z
#> 2 163007608 2025-05-19T14:31:37.000Z
#>                                  intended_use
#> 1 Primary routine opt-in example POD5 object.
#> 2                   Fail-state examples only.
```

Opt in explicitly before downloading from S3, and keep cached files
outside the repository:

``` r

pass_example <- ont_zymo_pod5_example_objects(role = "pass")

if (identical(tolower(Sys.getenv("FLOUNDER_RUN_NETWORK_TESTS")), "true")) {
  fetched <- ont_open_data_fetch(
    key = pass_example$key,
    cache_dir = tools::R_user_dir("floundeR", which = "cache")
  )

  pod5_file_info(fetched$cache_path)
  pod5_verify(fetched$cache_path)
}
```

## Local POD5 Inventory

The curated R API is intentionally small: discover POD5-containing
folders, verify candidate files, inspect file or folder metadata, build
manifests, compare collections, and plan read-only demonstration
subdivisions. These calls go through in-process Rust and return R data
frames or lists.

The tiny files below are POD5-like signature fixtures for documentation.
They exercise discovery, manifest, and envelope checks without embedding
real raw signal data in the package.

``` r

if (rust_available) {
  signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
  run_dir <- file.path(tempdir(), "flounder-pod5-qc-example")
  unlink(run_dir, recursive = TRUE)
  dir.create(file.path(run_dir, "pass"), recursive = TRUE)
  dir.create(file.path(run_dir, "fail"), recursive = TRUE)

  writeBin(
    c(signature, as.raw(rep(0, 16)), signature),
    file.path(run_dir, "pass", "example-pass.pod5")
  )
  writeBin(
    c(signature, as.raw(rep(1, 16)), signature),
    file.path(run_dir, "fail", "example-fail.pod5")
  )

  pod5_find(run_dir)
  pod5_folder_info(run_dir)
  pod5_manifest(run_dir)
} else {
  data.frame(
    status = "not_checked",
    reason = "Compiled Rust support is unavailable in this R installation.",
    stringsAsFactors = FALSE
  )
}
#>   schema_version
#> 1              1
#> 2              1
#>                                                                                  source
#> 1 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//Rtmp7TvK0e/flounder-pod5-qc-example
#> 2 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//Rtmp7TvK0e/flounder-pod5-qc-example
#>            relative_path
#> 1 fail/example-fail.pod5
#> 2 pass/example-pass.pod5
#>                                                                                                           path
#> 1 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//Rtmp7TvK0e/flounder-pod5-qc-example/fail/example-fail.pod5
#> 2 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//Rtmp7TvK0e/flounder-pod5-qc-example/pass/example-pass.pod5
#>   size_bytes verification_status verification_failed_checks
#> 1         32          incomplete                          0
#> 2         32          incomplete                          0
```

File-level inspection separates fast envelope evidence from unavailable
deep POD5 metadata. That prevents a lightweight tutorial fixture from
being overstated as fully parsed raw-signal evidence.

``` r

if (rust_available) {
  example_file <- file.path(run_dir, "pass", "example-pass.pod5")

  pod5_verify(example_file)
  pod5_file_info(example_file)
} else {
  data.frame(
    status = "not_checked",
    reason = "Compiled Rust support is unavailable in this R installation.",
    stringsAsFactors = FALSE
  )
}
#>                                                                                                           path
#> 1 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//Rtmp7TvK0e/flounder-pod5-qc-example/pass/example-pass.pod5
#>   size_bytes flow_cell_id sequencing_kit read_count acquisition_start_utc
#> 1         32         <NA>           <NA>         NA                  <NA>
#>   duration_seconds pod5_version integrity_status
#> 1               NA         <NA>      unavailable
#>                                                             integrity_reason
#> 1 POD5 parser backend not configured; only filesystem metadata was inspected
```

## Demonstration Subdivision Plans

Large POD5 files are not suitable for routine vignette builds. When a
real POD5 object has been cached outside the repository, floundeR can
ask the Rust-backed POD5 layer for a read-only split plan. This plans
demonstration subsets without writing new POD5 files.

``` r

if (rust_available) {
  pod5_subdivide_plan(run_dir, strategy = "file-count", files_per_chunk = 1)
} else {
  data.frame(
    status = "not_checked",
    reason = "Compiled Rust support is unavailable in this R installation.",
    stringsAsFactors = FALSE
  )
}
#>   schema_version
#> 1              1
#> 2              1
#>                                                                                  source
#> 1 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//Rtmp7TvK0e/flounder-pod5-qc-example
#> 2 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//Rtmp7TvK0e/flounder-pod5-qc-example
#>     strategy              target chunk_index chunk_label file_count total_bytes
#> 1 file-count 1 file(s) per chunk           1  chunk-0001          1          32
#> 2 file-count 1 file(s) per chunk           2  chunk-0002          1          32
#>   read_count         relative_paths warnings
#> 1         NA fail/example-fail.pod5     <NA>
#> 2         NA pass/example-pass.pod5     <NA>
```

Actual subset or split file creation remains an explicit derived-data
workflow: record the source URI, checksum, command/runtime provenance,
and generated manifest, then keep the derived POD5 artifact outside the
package repository. The maintained preparation contract is documented in
`DERIVED_POD5_DEMOS.md`.

## QC Handoff Shape

POD5 evidence is normally passed forward as manifest, verification, and
open-data provenance tables. Those tables can be attached to Synoptikon
handoff payloads or converted into Grammateus report elements when the
optional private runtime is installed.

``` r

pass_example <- ont_zymo_pod5_example_objects(role = "pass")

pod5_provenance <- data.frame(
  schema_version = "flounder.pod5_open_data_example.v1",
  status = "not_downloaded",
  bucket = pass_example$bucket,
  region = pass_example$region,
  key = pass_example$key,
  s3_uri = pass_example$s3_uri,
  size_bytes = pass_example$size,
  intended_use = pass_example$intended_use,
  stringsAsFactors = FALSE
)

pod5_provenance
#>                       schema_version         status        bucket    region
#> 1 flounder.pod5_open_data_example.v1 not_downloaded ont-open-data eu-west-1
#>                                                                              key
#> 1 zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_pass_279c9095_68316534_8289.pod5
#>                                                                                              s3_uri
#> 1 s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_pass_279c9095_68316534_8289.pod5
#>   size_bytes                                intended_use
#> 1   47077200 Primary routine opt-in example POD5 object.
```

Use the companion Synoptikon handoff vignette for the full JSON payload
contract.
