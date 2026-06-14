# Fetch one ONT open-data object into an explicit cache directory

`ont_open_data_fetch()` downloads exactly one selected public ONT
open-data object. It never expands a prefix into multiple downloads.
Call
[`ont_open_data_list()`](https://sagrudd.github.io/floundeR/reference/ont_open_data_list.md)
first when object discovery is required, then pass a single object key
to this function.

## Usage

``` r
ont_open_data_fetch(
  key,
  cache_dir = tools::R_user_dir("floundeR", which = "cache"),
  bucket = "ont-open-data",
  region = "eu-west-1",
  file_name = NULL,
  overwrite = FALSE,
  anonymous = TRUE
)
```

## Arguments

- key:

  Character scalar. Full S3 object key to download.

- cache_dir:

  Character scalar. Directory where the object should be cached.
  Defaults to the user cache for `floundeR`.

- bucket:

  Character scalar. Public S3 bucket name.

- region:

  Character scalar. AWS region for the bucket.

- file_name:

  Optional character scalar. Local file name inside `cache_dir`.
  Defaults to `basename(key)`.

- overwrite:

  Logical scalar. If `FALSE` and the local file already exists, return
  metadata without downloading again.

- anonymous:

  Logical scalar. If `TRUE`, temporarily sets `AWS_NO_SIGN_REQUEST=true`
  while listing and downloading.

## Value

A one-row data frame with `bucket`, `key`, `size`, `last_modified_utc`,
`etag`, `storage_class`, `cache_path`, and `downloaded`.

## Examples

``` r
if (FALSE) { # \dontrun{
if (identical(tolower(Sys.getenv("FLOUNDER_RUN_NETWORK_TESTS")), "true")) {
  ont_open_data_fetch(
    key = paste0(
      flounder_ont_zymo_pod5_prefix(),
      "PAU85136_pass_279c9095_68316534_8289.pod5"
    ),
    cache_dir = file.path(tempdir(), "flounder-ont-open-data")
  )
}
} # }
```
