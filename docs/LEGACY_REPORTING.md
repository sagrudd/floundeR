# Legacy RMarkdown Reporting Boundary

`floundeR` no longer treats RMarkdown as the primary report-rendering
stack for nanopore QC. Governed technical reports are expected to flow
through Grammateus semantic report contracts, optional private
Grammateus runtime assets, and Mnemosyne Biosciences templates/themes.

RMarkdown support is retained only for transitional package
documentation while the reboot reaches full Grammateus report coverage.

## Allowed RMarkdown Uses

- Package vignettes that teach open-source API usage.
- Transitional examples that explain how existing R workflows map into
  Grammateus contracts.
- Legacy migrated tutorials kept outside the source-package build.
- Non-governed developer notes during the reboot, provided they are not
  presented as production report output.

## Disallowed RMarkdown Uses

- Production QC report rendering.
- Mnemosyne Biosciences branded technical reports.
- Synoptikon trusted-report lifecycle artifacts.
- New APIs that return RMarkdown documents as the normal report format.
- Parallel CSS, template, signing, provenance, or approval models that
  bypass Grammateus.

## Current Repository State

- `vignettes/synoptikon-handoff.Rmd` remains an R package vignette. It
  is a network-free handoff example, not the governed report renderer.
- `legacy-vignettes/` contains historical RMarkdown material and is
  excluded from source-package builds.
- `docs/` contains release-oriented pkgdown output and is not the source
  of truth for current report behavior.
- `rmarkdown` and `knitr` remain in `Suggests`/`VignetteBuilder` while
  package vignettes are still built through the standard R toolchain.

## Migration Rule

When a reporting example becomes part of the normal QC report surface,
express it as:

1.  stable QC evidence tables;
2.  [`grammateus_report_element()`](https://sagrudd.github.io/floundeR/reference/grammateus_report_element.md)
    or
    [`grammateus_qc_report_elements()`](https://sagrudd.github.io/floundeR/reference/grammateus_report_element.md)
    objects;
3.  governed `flounder_grammateus_figure` or
    `flounder_grammateus_plot_spec` objects;
4.  [`qc_report()`](https://sagrudd.github.io/floundeR/reference/qc_report.md)
    contract and manifest output;
5.  optional HTML/PDF rendering through a validated Grammateus runtime.

RMarkdown may describe that workflow, but it must not own the workflow.
