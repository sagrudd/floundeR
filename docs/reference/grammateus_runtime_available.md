# Discover and validate a prebuilt Grammateus runtime

These helpers inspect optional private Grammateus runtime bundles
without making Grammateus a required dependency of the public floundeR
package. Discovery checks `GRAMMATEUS_HOME`, then
`options(floundeR.grammateus_home = "...")`, then the floundeR user
cache. No helper downloads private assets during package load or
ordinary QC work.

## Usage

``` r
grammateus_runtime_available(runtime_root = NULL)

grammateus_runtime_version(runtime_root = NULL)

grammateus_runtime_manifest(runtime_root = NULL)

grammateus_runtime_validate(
  runtime_root = NULL,
  required_capabilities = c("render_report_html", "render_report_pdf")
)

grammateus_runtime_install(
  source_runtime_root,
  cache_root = tools::R_user_dir("floundeR", "cache"),
  overwrite = FALSE
)
```

## Arguments

- runtime_root:

  Optional explicit runtime root. When omitted, the standard discovery
  order is used.

- required_capabilities:

  Character vector of manifest capabilities that must be `TRUE`.

- source_runtime_root:

  Existing extracted runtime root to copy into the floundeR-managed
  cache.

- cache_root:

  Optional cache root. Defaults to
  `tools::R_user_dir("floundeR", "cache")`.

- overwrite:

  Whether `grammateus_runtime_install()` may replace an existing cached
  runtime directory.

## Value

`grammateus_runtime_available()` returns a logical scalar.
`grammateus_runtime_version()` returns a character scalar or
`NA_character_`. `grammateus_runtime_manifest()` returns the parsed
manifest list or `NULL`. `grammateus_runtime_validate()` returns a
validation list with schema version, availability, validity, runtime
path, version, platform, capabilities, artifact counts, aggregate
artifact hash, and failure records. `grammateus_runtime_install()`
returns the validation list for the installed runtime.
