# floundeR Dependency Assessment

This file records the current R dependency posture for the revival
build. It is generated from the declared `DESCRIPTION` dependency
surface and should be refreshed when dependencies are added, removed, or
moved between `Imports` and `Suggests`.

## Container Baseline

The development Docker image builds from `rocker/r-ver:4.6.0` and
installs the declared CRAN and Bioconductor dependency set through
`scripts/bootstrap-r-dependencies.R --install`.

Verified on 2026-06-13:

- R: 4.6.0
- Bioconductor: 3.23
- Image tag: `flounder-dev:latest`
- Build result: passed
- Notable native dependency coverage: `Rsamtools` and `ShortRead` built
  successfully from source with the system libraries declared in
  `Dockerfile`.

The image is intentionally a development/check image. It includes Cargo
and rustc so source package checks can build the embedded Rust
extension. It does not include private Grammateus runtime assets,
downloaded POD5 files, or derived large example data.

## Version Snapshot

The following versions were reported by `scripts/audit-r-dependencies.R`
inside the clean container after dependency installation.

| Package | Declared In | Channel | Version | Support Posture |
|----|----|----|----|----|
| BiocGenerics | Imports | Bioconductor | 0.58.1 | Manage through `BiocManager` for the active R/Bioconductor release. |
| Biostrings | Imports | Bioconductor | 2.80.1 | Supported Bioconductor sequence infrastructure; retained for FASTA/FASTQ and sequence evidence while modern Rust-backed paths are added. |
| GenomeInfoDb | Imports | Bioconductor | 1.48.0 | Supported Bioconductor genome metadata infrastructure; review necessity as core QC contracts narrow. |
| GenomicRanges | Imports | Bioconductor | 1.64.0 | Supported Bioconductor ranges infrastructure; review necessity as BAM/POD5 Rust bindings mature. |
| IRanges | Imports | Bioconductor | 2.46.0 | Supported Bioconductor ranges infrastructure; dependency burden should remain justified by active code. |
| Rsamtools | Imports | Bioconductor | 2.28.0 | Supported but native-heavy; keep only until curated Bamana-backed BAM QC can replace or narrow its role. |
| ShortRead | Imports | Bioconductor | 1.70.0 | Supported but native-heavy and legacy-oriented; retain only for existing FASTQ/FASTA paths until lighter QC wrappers exist. |
| aws.s3 | Imports | CRAN | 0.3.22 | Retained for anonymous ONT open-data access; network tests remain opt-in. |
| cli | Imports | CRAN | 3.6.6 | Supported CRAN diagnostics dependency. |
| dplyr | Imports | CRAN | 1.2.1 | Supported tidy data dependency; should remain if public QC outputs stay tidy. |
| emojifont | Imports | CRAN | 0.6.0 | Legacy presentation dependency; candidate for removal unless active report styling still requires it. |
| ggplot2 | Imports | CRAN | 4.0.3 | Supported plotting dependency; central for R plot handoff into Grammateus. |
| ggsci | Imports | CRAN | 5.0.0 | Palette dependency; review against Mnemosyne/Grammateus branding before retaining. |
| jsonlite | Imports | CRAN | 2.0.0 | Supported JSON dependency; required for Synoptikon QC payload writing and schema-facing tests. |
| knitr | Suggests | CRAN | 1.51 | Standard R package vignette builder; not used as the production QC report renderer. |
| lubridate | Imports | CRAN | 1.9.5 | Supported timestamp handling dependency. |
| magick | Imports | CRAN | 2.9.1 | Native image dependency; review once Grammateus asset rendering owns report images. |
| magrittr | Imports | CRAN | 2.0.5 | Pipe dependency used by legacy API style; review after modernisation. |
| purrr | Imports | CRAN | 1.2.2 | Supported tidy iteration dependency. |
| R6 | Imports | CRAN | 2.6.1 | Legacy object model dependency; retain while existing classes remain exported. |
| RColorBrewer | Imports | CRAN | 1.1-3 | Palette dependency; review against current branding. |
| readr | Imports | CRAN | 2.2.0 | Supported table parser; central to sequencing-summary modernisation. |
| reshape2 | Imports | CRAN | 1.4.5 | Legacy reshape dependency; candidate for replacement with tidyr/dplyr patterns. |
| rlang | Imports | CRAN | 1.2.0 | Supported tidy error/programming dependency. |
| rmarkdown | Suggests | CRAN | 2.31 | Standard R package vignette engine only; Grammateus remains the strategic HTML/PDF QC report renderer. |
| stringdist | Imports | CRAN | 0.9.17 | Supported fuzzy matching dependency; likely relevant to kit/adapter evidence. |
| stringr | Imports | CRAN | 1.6.0 | Supported string handling dependency. |
| testthat | Suggests | CRAN | 3.3.2 | Supported test framework dependency. |
| tibble | Imports | CRAN | 3.3.1 | Supported tidy return-value dependency. |
| tidyr | Imports | CRAN | 1.3.2 | Supported tidy reshaping dependency. |
| tidyselect | Imports | CRAN | 1.2.1 | Supported tidy selection dependency. |

## Current Check Findings

After the container dependency build,
`R CMD check --no-manual --no-build-vignettes .` progressed beyond
dependency validation. The remaining issues are package-quality tasks,
not missing dependency blockers:

- `Blast` examples/tests fail because
  `flnDr("drosophila_uniref100.blastx.gz")` resolves to `NA` under
  check.
- The POD5 manifest fixture has at least one `sha256` value that is not
  64 characters.
- Legacy vignettes fail during direct source-tree check execution.
- Direct source-tree checks report expected notes/warnings for hidden
  development files; release checks should be run on an `R CMD build`
  tarball.
- Several declared imports are currently unused by static analysis and
  should be removed, moved, or justified during API modernisation.

## Commands

Build the dependency image:

``` sh
docker build -t flounder-dev .
```

Audit dependencies in the image:

``` sh
docker run --rm flounder-dev Rscript scripts/audit-r-dependencies.R
```

Run the source-tree package check in the image:

``` sh
docker run --rm flounder-dev R CMD check --no-manual --no-build-vignettes .
```
