# Verify a local POD5 candidate file

`pod5_verify()` checks a single local file through the in-process Rust
extension and the curated `pod5-tools` library API. The current backend
performs fast extension and fixed POD5 signature checks, then reports
deeper layout/schema checks as `not_checked` until the full parser
backend is connected. A signature-valid file therefore returns
`overall_status = "incomplete"` rather than overstating full POD5
validity.

## Usage

``` r
pod5_verify(path)
```

## Arguments

- path:

  Character scalar. Candidate POD5 file to verify.

## Value

A data frame with one row per verification check and the columns `path`,
`size_bytes`, `overall_status`, `check`, `category`, `status`, and
`detail`.

## Examples

``` r
signature <- as.raw(c(0x8b, 0x50, 0x4f, 0x44, 0x0d, 0x0a, 0x1a, 0x0a))
path <- file.path(tempdir(), "flounder-signature-example.pod5")
writeBin(c(signature, as.raw(rep(0, 16)), signature), path)
pod5_verify(path)
#>                                                                                           path
#> 1 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//RtmpfBxQW1/flounder-signature-example.pod5
#> 2 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//RtmpfBxQW1/flounder-signature-example.pod5
#> 3 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//RtmpfBxQW1/flounder-signature-example.pod5
#> 4 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//RtmpfBxQW1/flounder-signature-example.pod5
#> 5 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//RtmpfBxQW1/flounder-signature-example.pod5
#> 6 /var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T//RtmpfBxQW1/flounder-signature-example.pod5
#>   size_bytes overall_status                check  category      status
#> 1         32     incomplete            extension extension      passed
#> 2         32     incomplete    leading_signature signature      passed
#> 3         32     incomplete   trailing_signature signature      passed
#> 4         32     incomplete combined_file_layout    layout not_checked
#> 5         32     incomplete      required_tables    schema not_checked
#> 6         32     incomplete      schema_metadata    schema not_checked
#>                                                                                                                   detail
#> 1                                                                                                file extension is .pod5
#> 2                                                                           leading signature matches ONT POD5 signature
#> 3                                                                          trailing signature matches ONT POD5 signature
#> 4 combined-file layout, section markers, footer magic, footer length, and padding checks require the POD5 parser backend
#> 5                                      Reads, Signal, and Run Info table presence checks require the POD5 parser backend
#> 6                  POD5 version, writer software, and file identifier consistency checks require the POD5 parser backend
```
