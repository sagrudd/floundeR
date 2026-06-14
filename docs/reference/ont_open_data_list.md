# List objects in an ONT open-data S3 prefix

`ont_open_data_list()` performs an explicit anonymous listing of a
public ONT open-data S3 prefix. The default prefix is the selected Zymo
fecal POD5 source used by floundeR integration examples. The helper
lists metadata only; it does not download POD5 files.

## Usage

``` r
ont_open_data_list(
  prefix = flounder_ont_zymo_pod5_prefix(),
  bucket = "ont-open-data",
  region = "eu-west-1",
  max = 1000L,
  anonymous = TRUE
)

flounder_ont_zymo_pod5_prefix()
```

## Arguments

- prefix:

  Character scalar. S3 object prefix to list.

- bucket:

  Character scalar. Public S3 bucket name.

- region:

  Character scalar. AWS region for the bucket.

- max:

  Integer scalar. Maximum number of objects to return.

- anonymous:

  Logical scalar. If `TRUE`, temporarily sets `AWS_NO_SIGN_REQUEST=true`
  while listing.

## Value

A data frame with one row per object and the columns `bucket`, `key`,
`size`, `last_modified_utc`, `etag`, and `storage_class`.

## Examples

``` r
if (FALSE) { # \dontrun{
if (identical(tolower(Sys.getenv("FLOUNDER_RUN_NETWORK_TESTS")), "true")) {
  ont_open_data_list(max = 5)
}
} # }
```
