# Build and render Grammateus semantic plot specifications

`grammateus_plot_spec()` creates a backend-neutral Grammateus
`ReportPlot`-shaped object from tidy R data. `grammateus_render_plot()`
executes the controlled R/ggplot2 backend for that spec, writing
deterministic run inputs, a generated script, backend session metadata,
and requested PNG/SVG artifacts. These helpers keep controlled plot
generation available from the open-source package without requiring
private Grammateus runtime assets.

## Usage

``` r
grammateus_plot_spec(
  plot_id,
  plot_type = c("line", "bar", "stacked_bar", "scatter", "pca"),
  data,
  mappings,
  axes,
  caption,
  theme = "mnemosyne_qc",
  palette = "viridis",
  output = list(width_mm = 160, height_mm = 100, dpi = 300, formats = c("svg",
    "png")),
  bar_value_semantics = NULL,
  source_hash = NULL,
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL
)

grammateus_render_plot(
  plot_spec,
  execution = c("local_rscript", "docker_container"),
  run_root,
  rscript_path = "Rscript",
  docker_path = "docker",
  container_image = NULL,
  container_rscript = "Rscript",
  required_packages = c("ggplot2", "jsonlite", "scales", "svglite", "viridisLite"),
  overwrite_run_dir = TRUE
)
```

## Arguments

- plot_id:

  Stable lower-snake-case plot identifier. Must start with `plot_`.

- plot_type:

  Semantic plot family: `line`, `bar`, `stacked_bar`, `scatter`, or
  `pca`.

- data:

  Tidy data frame for inline plot data, or a Grammateus data reference
  list with `kind = "reference"`, `path`, `source_hash`, and
  `row_count`.

- mappings:

  List with required `x` and `y` fields plus optional `color`, `fill`,
  `group`, and `label` fields.

- axes:

  List with `x` and `y` axis definitions. Each axis requires `label` and
  may include `unit`, `scale`, and `variance_explained`.

- caption:

  Scientific plot caption.

- theme:

  Named Grammateus theme profile.

- palette:

  Named Grammateus palette policy.

- output:

  Output definition with `width_mm`, `height_mm`, `dpi`, and `formats`.

- bar_value_semantics:

  Required for `bar` and `stacked_bar` plots. One of `counts`,
  `percentages`, `rates`, or `normalized_measurements`.

- source_hash:

  Optional SHA-256 hash for the source data. Plain 64-character hex
  digests are accepted and normalised to `sha256:<hex>`.

- produced_by:

  Producing software or service name.

- producer_version:

  Producing software version. When `NULL`, the installed floundeR
  version is used.

- produced_at_utc:

  Production timestamp.

- run_id:

  Optional upstream run, workflow, or analysis identifier.

- plot_spec:

  A `flounder_grammateus_plot_spec` object.

- execution:

  Plot backend execution mode. `local_rscript` runs an `Rscript`
  executable directly. `docker_container` runs `Rscript` inside a
  Docker-compatible image with the run directory mounted at `/work`.

- run_root:

  Directory where deterministic per-plot run directories are created.

- rscript_path:

  Path to the local `Rscript` executable.

- docker_path:

  Path to the Docker-compatible CLI executable.

- container_image:

  Container image used when `execution` is `docker_container`.

- container_rscript:

  Rscript command inside the container.

- required_packages:

  R packages required by the generated script.

- overwrite_run_dir:

  Whether an existing deterministic run directory may be replaced.

## Value

`grammateus_plot_spec()` returns a list with class
`flounder_grammateus_plot_spec`. The object follows the Grammateus
`ReportPlot` shape and contains plot identity, type, data, mappings,
axes, output settings, optional bar-value semantics, and provenance.

`grammateus_render_plot()` returns a list with class
`flounder_grammateus_plot_run` containing run paths, hashes, stdout,
stderr, backend session path, and artifact metadata. Artifacts include
`flounder_grammateus_figure` wrappers for the generated PNG/SVG files.
