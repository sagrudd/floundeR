# Describe the canonical ONT Zymo fecal POD5 dataset

`ont_zymo_pod5_dataset()` returns stable metadata for the public ONT
Zymo fecal POD5 prefix used by floundeR's opt-in real-data
demonstrations. It does not contact S3.

## Usage

``` r
ont_zymo_pod5_dataset()
```

## Value

A one-row data frame with dataset name, bucket, region, prefix, S3 URI,
observed object counts, observed total bytes, and verification date.

## Examples

``` r
ont_zymo_pod5_dataset()
#>                           dataset_name        bucket    region
#> 1 ont_zymo_fecal_2025_05_pau85136_pod5 ont-open-data eu-west-1
#>                                  prefix
#> 1 zymo_fecal_2025.05/raw/PAU85136/pod5/
#>                                                     s3_uri total_pod5_objects
#> 1 s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/               9107
#>   pass_pod5_objects fail_pod5_objects  total_bytes         verified_utc
#> 1              8290               817 2.676974e+12 2026-06-13T00:00:00Z
```
