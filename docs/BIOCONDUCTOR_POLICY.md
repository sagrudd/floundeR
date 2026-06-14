# Bioconductor Policy Posture

This document records the Bioconductor-facing policy check for floundeR
before any future submission. It was verified against official
Bioconductor package guidance on 2026-06-14. Bioconductor documentation
cited here was last built on 2026-06-09.

floundeR is not being submitted to Bioconductor in this slice. This is a
readiness posture for later release planning.

## Sources Checked

- Bioconductor package `DESCRIPTION` guidance:
  <https://contributions.bioconductor.org/description.html>
- Bioconductor README guidance:
  <https://contributions.bioconductor.org/readme.html>
- Bioconductor package data guidance:
  <https://contributions.bioconductor.org/data.html>
- Bioconductor R-code guidance, including additional files,
  dependencies, web queries, and file caching:
  <https://contributions.bioconductor.org/r-code.html>
- Bioconductor unit test guidance:
  <https://contributions.bioconductor.org/tests.html>
- Bioconductor third-party/native-code guidance:
  <https://contributions.bioconductor.org/other-than-Rcode.html>
- Bioconductor build/check/BiocCheck guidance:
  <https://contributions.bioconductor.org/build-check-bioccheck.html>
- Bioconductor web-resource guidance:
  <https://contributions.bioconductor.org/querying-web-resources.html>

## Policy Conclusions

### Private Grammateus Runtime Assets

Bioconductor-facing floundeR must not require private Grammateus source
or private Grammateus runtime downloads for installation, package load,
examples, vignettes, tests, `R CMD check`, or `BiocCheck`.

The public package may document that governed HTML/PDF rendering is
enhanced by an optional private runtime, but that runtime must remain
outside the Bioconductor package and outside default checks. Core QC,
report-contract creation, plot/figure preparation, and Synoptikon
payload generation must continue to work without Grammateus assets.

Required floundeR posture:

- Keep Grammateus runtime support optional and explicit.
- Keep runtime downloads out of package code, examples, README-evaluated
  code, vignettes, and tests.
- Keep private runtime integration tests opt-in through environment
  variables and an explicit `GRAMMATEUS_HOME`.
- Keep public tests asserting that core QC works when Grammateus is
  absent.
- Do not include private Grammateus source, runtime archives, extracted
  runtime assets, credentials, or private release URLs in the submitted
  package.

### External System Requirements

Bioconductor permits external software requirements when they are
documented in `SystemRequirements`, but package code must not install
system dependencies for users. Non-trivial installation steps should be
described in a top-level `INSTALL` file and user-facing documentation.

Required floundeR posture:

- Keep Rust/Cargo requirements in `SystemRequirements`.
- Add or maintain a top-level `INSTALL` file before submission, covering
  Linux, macOS, and Windows source-install requirements.
- Do not install Rust, Cargo, system libraries, private runtimes, or
  external applications from R code, tests, examples, vignettes, or
  README-evaluated chunks.
- Prefer open, trusted, easy-to-install system requirements.
- Treat Grammateus runtime assets as optional external artifacts, not
  required system software.

### R Package Dependencies

Bioconductor package dependencies must be available through Bioconductor
or CRAN. The `Remotes:` field is not supported for Bioconductor builds.

Required floundeR posture:

- Do not depend on GitHub-only R packages for any Bioconductor-facing
  build.
- Keep optional developer-only tools out of mandatory dependency fields.
- Ensure all `Imports`, `Depends`, and `Suggests` packages are available
  from CRAN or Bioconductor before submission.
- Do not encode private Grammateus source as an R package dependency.

### Downloads, Web Resources, And Caching

Bioconductor allows web-resource access when it is robust, bounded,
clear on failure, and does not make checks slow or unreliable.
Downloaded files should use Bioconductor-friendly caching, preferably
`BiocFileCache`, or standard cache locations such as
[`tools::R_user_dir()`](https://rdrr.io/r/tools/userdir.html) when
appropriate. Package code must not download or write files to a user’s
home directory, working directory, or installed package directory by
default.

Required floundeR posture:

- Keep ONT S3 examples and POD5 downloads opt-in.
- Keep network tests skipped by default.
- Keep default package checks runnable without network access.
- Cache external files outside the repository and outside the installed
  package directory.
- Prefer `BiocFileCache` for Bioconductor-facing download helpers, or
  document why `tools::R_user_dir("floundeR", "cache")` is used for a
  specific cache.
- Fail quickly and clearly when web resources are unavailable.
- Keep default examples and vignettes under reasonable check-time
  limits.

### Large Data And POD5 Examples

Bioconductor does not support git-lfs for package data, and large
datasets should be routed through ExperimentHub or AnnotationHub where
appropriate. Smaller raw data fixtures may live under `inst/extdata`
when documented with source/provenance scripts.

Required floundeR posture:

- Do not commit downloaded ONT POD5 files or large derived POD5
  artifacts.
- Keep the ONT Zymo fecal POD5 source as opt-in real-data documentation,
  not a default package fixture.
- Before submission, decide whether any real POD5 example should be
  represented by a tiny documented fixture, an
  ExperimentHub/AnnotationHub-backed resource, or a non-running opt-in
  example.
- Keep source object keys, sizes, timestamps, checksums, and derivation
  provenance for any example data workflow.

### Tests And Checks

Bioconductor expects packages to build, check, and pass `BiocCheck`
without errors. Unit tests are part of the standard build/check process.

Required floundeR posture:

- Run `R CMD build`, `R CMD check`, and `BiocCheck` against the current
  Bioconductor devel environment before submission.
- Keep tests offline by default.
- Keep private Grammateus tests behind explicit opt-in flags and runtime
  paths.
- Keep tests for missing optional runtime behavior in the public test
  suite.
- Reassess long-running or network-sensitive tests before submission.

## Current floundeR Readiness Assessment

Current posture is directionally compatible with a future Bioconductor
track because:

- Grammateus remains optional and private.
- Core QC APIs are designed to work without Grammateus runtime assets.
- ONT POD5 real-data usage is opt-in and cache-oriented.
- GitHub Actions are currently disabled, but disabled templates separate
  public runtime-free checks from private credentialed runtime checks.

Known submission blockers or follow-up work:

- Add a top-level `INSTALL` file for source-install system requirements.
- Add `biocViews` before any formal submission.
- Confirm all package dependencies are CRAN/Bioconductor-available and
  remove any GitHub-only dependency assumptions.
- Prefer or support `BiocFileCache` for Bioconductor-facing external
  file caching.
- Run `BiocCheck` in a current Bioconductor devel environment.
- Resolve current local dependency gaps that prevent full `R CMD check`
  in this host library.
- Keep private Grammateus runtime assets, credentials, URLs, and
  workflows out of the public package build path.
