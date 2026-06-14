# Return selected ONT Zymo fecal POD5 example objects

`ont_zymo_pod5_example_objects()` returns the fixed pass/fail POD5
object choices used for floundeR examples and opt-in integration
workflows. The pass object is the routine example. The fail object
should be used only where fail-state QC evidence is required. This
helper does not contact S3.

## Usage

``` r
ont_zymo_pod5_example_objects(role = c("all", "pass", "fail"))
```

## Arguments

- role:

  Character scalar. Return all selected examples, only the pass example,
  or only the fail-state example.

## Value

A data frame with `role`, `state`, `file_name`, `bucket`, `region`,
`key`, `s3_uri`, `size`, `last_modified_utc`, and `intended_use`.

## Examples

``` r
ont_zymo_pod5_example_objects()
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
ont_zymo_pod5_example_objects(role = "pass")$key
#> [1] "zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_pass_279c9095_68316534_8289.pod5"
```
