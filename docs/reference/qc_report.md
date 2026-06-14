# Assemble a Grammateus-backed QC report contract

`qc_report()` is the high-level floundeR report assembly API. It
combines prepared Grammateus semantic report elements, optional governed
figures, and a Mnemosyne Biosciences theme descriptor into a stable
report contract and manifest. The public package always writes the
contract and manifest without requiring private Grammateus assets.
HTML/PDF rendering is attempted only when requested and an authorized
Grammateus runtime is available; otherwise the manifest records an
explicit render status.

## Usage

``` r
qc_report(
  elements,
  figures = NULL,
  output_dir,
  output = c("html", "pdf"),
  render = c("if_available", "never", "require"),
  theme = grammateus_mnemosyne_theme(),
  report_id = "report_nanopore_qc",
  title = "Nanopore sequencing QC report",
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  overwrite = TRUE
)
```

## Arguments

- elements:

  A single Grammateus report element, a
  `flounder_grammateus_report_element_bundle`, or a named list of report
  elements.

- figures:

  Optional governed figure, figure bundle, or list of governed figures
  prepared by
  [`grammateus_figure_from_file()`](https://sagrudd.github.io/floundeR/reference/grammateus_figure_from_file.md)
  or
  [`grammateus_figure_from_ggplot()`](https://sagrudd.github.io/floundeR/reference/grammateus_figure_from_file.md).

- output_dir:

  Directory where the report contract and manifest are written.

- output:

  Requested rendered formats. Use any of `html` and `pdf`.

- render:

  Render policy: `if_available` records unavailable render outputs when
  the private runtime is absent, `never` writes only the contract and
  manifest, and `require` errors unless rendering is available.

- theme:

  A `flounder_grammateus_theme` object.

- report_id:

  Stable lower-snake-case report identifier. Must start with `report_`.

- title:

  Human-readable report title.

- produced_by:

  Producing software or service name.

- producer_version:

  Producing software version.

- produced_at_utc:

  Production timestamp.

- run_id:

  Optional upstream run, workflow, or analysis identifier.

- overwrite:

  Whether existing report artifact paths may be replaced.

## Value

A list with class `flounder_qc_report` containing the themed report
contract, manifest, contract path, manifest path, requested output
status, and provenance.
