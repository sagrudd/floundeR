# Build Mnemosyne Grammateus report theme metadata

`grammateus_mnemosyne_theme()` returns a runtime-free descriptor for the
Mnemosyne Biosciences Grammateus template and theme expected by floundeR
QC reports. It does not bundle private Grammateus assets or inline CSS;
it records the template, theme, brand, style policy, and provenance that
an authorized Grammateus runtime can resolve during rendering.

## Usage

``` r
grammateus_mnemosyne_theme(
  profile = c("technical_qc", "appendix"),
  palette = c("mnemosyne_qc"),
  runtime_theme_id = "mnemosyne_biosciences_qc_v1",
  template_id = "mnemosyne_technical_qc_v1",
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL
)

grammateus_apply_theme(
  elements,
  theme = grammateus_mnemosyne_theme(),
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL
)
```

## Arguments

- profile:

  Mnemosyne report profile.

- palette:

  Mnemosyne palette policy.

- runtime_theme_id:

  Grammateus runtime theme identifier.

- template_id:

  Grammateus runtime template identifier.

- produced_by:

  Producing software or service name.

- producer_version:

  Producing software version.

- produced_at_utc:

  Production timestamp.

- run_id:

  Optional upstream run, workflow, or analysis identifier.

- elements:

  A single Grammateus report element, a
  `flounder_grammateus_report_element_bundle`, or a named list of report
  elements.

- theme:

  A `flounder_grammateus_theme` object.

## Value

`grammateus_mnemosyne_theme()` returns a list with class
`flounder_grammateus_theme`. `grammateus_apply_theme()` returns a list
with class `flounder_grammateus_themed_report`.

## Details

`grammateus_apply_theme()` wraps prepared Grammateus semantic elements
with a theme descriptor so report assembly can pass one coherent
themed-report contract to the rendering layer.
