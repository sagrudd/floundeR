# AGENTS.md

This repository contains `floundeR`, an R package for tidy analysis and QC of
Oxford Nanopore sequencing data. The reboot direction is to make `floundeR` a
contemporary R-native QC package with POD5-aware raw-data inspection backed by
in-process Rust code. The toolbox is not only a POD5 project; BAM/BGZF/FASTQ
processing remains a first-class QC surface.

## Project Direction

- Retire FAST5 capabilities and remove the `rhdf5` dependency.
- Do not implement normal POD5 support by shelling out to a Rust CLI.
- Use Rust inside R for POD5-heavy operations, preferably through `extendr` or
  an equivalent compiled R extension boundary.
- Reuse and extend `../pod5-tools` as a Rust library where appropriate, but keep
  R-facing APIs native: R functions should return data frames, tibbles, lists,
  and R conditions.
- Reuse and extend `../bamana` as a Rust library for BAM, BGZF, FASTQ,
  sampling, ingest, summary, validation, and forensic QC surfaces.
- Do not mirror every `pod5-tools` or `bamana` capability in `floundeR`.
  Integrate only the subset that makes `floundeR` a best-in-breed nanopore QC
  and review toolbox.
- Move report rendering away from RMarkdown-first workflows and toward
  Grammateus semantic report contracts, rendered through Rust embedded inside R.
- Keep `floundeR` fully open-source. Grammateus remains private, so it must be
  consumed as an optional prebuilt runtime/content bundle rather than a required
  public source dependency.
- Maintain Mnemosyne Biosciences branding, styling, provenance, and standard
  technical-report structure in generated reports.
- Prioritise QC contracts needed by synoptikon via `../mnemosyne`.
- It is permitted to modify `../pod5-tools` and `../grammateus` when required
  to provide clean Rust library APIs for `floundeR`.
- It is permitted to modify `../bamana` when required to provide clean Rust
  library APIs for `floundeR`.
- Cross-repository changes must not break product charters, SRS material, ADRs,
  ARDs, or canonical contracts in `../mnemosyne-docs`.
- `pod5-tools` and `bamana` are open-source integration engines; Grammateus is
  private runtime/reporting infrastructure. Do not blur those distribution
  boundaries.
- Use ONT open-data POD5 examples for integration and documentation, especially
  `s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/`, but never commit
  downloaded POD5 files.
- Use `PAU85136_pass_279c9095_68316534_8289.pod5` as the primary pass example
  and `PAU85136_fail_279c9095_68316534_0.pod5` only for fail-state examples.
- It is acceptable to use `../pod5-tools` to subset or split selected POD5
  source files for demonstrations, provided provenance and checksums are
  recorded and large derived artifacts are not committed to the package repo.

## Development Rules

- Read the relevant R code, tests, and package metadata before editing.
- Keep changes scoped to the requested slice.
- Prefer existing package style unless changing it is part of an explicit
  modernisation task.
- Do not leave generated documentation, namespace files, or tests inconsistent
  with source changes.
- Do not remove user changes from the working tree unless explicitly asked.

## R Package Standards

- Keep `DESCRIPTION`, `NAMESPACE`, `R/`, `man/`, tests, vignettes, and README in
  agreement.
- Avoid `eval(parse())`; use structured R APIs.
- Public functions must have documented return shapes.
- Public QC outputs should be stable and schema-versioned where synoptikon may
  consume them.
- Tests should run without network access and should skip clearly for optional
  system requirements.
- Network-dependent examples and integration tests must be opt-in.
- Large external datasets must be downloaded into an explicit cache directory
  outside the repository by default.

## Rust-In-R Standards

- Rust code used by R must expose library functions, not rendered CLI output.
- R wrappers should convert Rust structs into idiomatic R objects.
- Rust errors should become predictable R conditions.
- Keep process spawning out of the main package API.
- If `../pod5-tools` needs changes to become a clean library dependency, make
  those changes deliberately and test both sides.
- If `../bamana` needs changes to become a clean BAM-processing library
  dependency, make those changes deliberately and test both sides.
- If `../grammateus` needs changes to become a clean report-rendering library
  dependency, make those changes deliberately and test both sides.
- Before changing `../pod5-tools` or `../grammateus`, inspect the relevant
  canonical documents in `../mnemosyne-docs`, especially product charters and
  current contract documents.
- Before changing `../bamana`, inspect `../bamana/CHARTER.md`,
  `../bamana/docs/project-charter.md`, and the relevant public JSON/schema and
  command documentation. If canonical bamana docs later appear in
  `../mnemosyne-docs`, inspect those too.
- If implementation requirements conflict with canonical docs, stop and update
  the plan rather than silently diverging.

## BAM/Bamana Standards

- Do not reimplement BAM parsing in R when Bamana can own the operation.
- Do not shell out to the `bamana` CLI for normal floundeR APIs.
- Keep the floundeR BAM API curated around QC, review, provenance, reporting,
  and synoptikon handoff rather than exposing Bamana wholesale.
- Preserve Bamana's distinctions between verification, validation, EOF checks,
  index checks, sort checks, mapping summaries, header-only mutation,
  record-level read-group annotation, collection-duplication inspection,
  collection-duplication remediation, and provenance/forensic inspection.
- R wrappers should expose stable data frames/lists and predictable R
  conditions while preserving Bamana response semantics.
- Transformation operations such as filtering, subsampling, sorting, merging,
  unmapping, reheadering, and annotation must require explicit output paths and
  provenance reporting.

## Reporting Standards

- Grammateus is the target reporting engine for technical QC reports.
- Grammateus support must be optional from the open-source package's point of
  view: core QC should install and run without private Grammateus source or
  private runtime assets.
- RMarkdown may be retained only as a legacy or transitional path while
  Grammateus coverage is incomplete.
- Reporting must support both directions of control: R can drop existing
  `ggplot2` plot artifacts into a Grammateus report, and Grammateus can call
  controlled R/ggplot2 generation from semantic plot specs.
- Reports must use semantic elements for tables, figures, plots, images,
  provenance, methods, limitations, and appendices.
- Tables, figures, and plots must have stable identifiers and captions.
- Image-like elements must include alt text for HTML output.
- Reports must preserve Mnemosyne Biosciences branding and styling through
  shared templates/themes rather than one-off CSS.
- Generated report artifacts should include manifests and provenance suitable
  for Synoptikon trusted-report lifecycle handling.
- Reports that include aligned-read evidence should include Bamana-derived BAM
  summary, validation, index, mapping, sorting, tag, and provenance sections as
  appropriate.
- The detailed R/Grammateus interface is documented in
  `REPORTING_INTERFACE.md`; keep it current when report APIs change.
- The distribution model for private Grammateus runtime assets is documented in
  `DISTRIBUTION.md`; keep it current when release or install behavior changes.

## Semantic Versioning

This package must follow semantic versioning.

- Patch version bumps are expected for bug fixes, dependency cleanup, test
  fixes, documentation corrections that affect installation/use, and small
  backwards-compatible additions.
- Minor version bumps are expected for substantive new features, new public
  functions, new QC schemas, POD5 Rust integration, or planned deprecations.
- Major version bumps are reserved for deliberate breaking API changes after
  the reboot reaches a stable release line.

Agents are expected to bump the patch or minor version after substantive
prompts, and to add a corresponding `NEWS.md` entry. If a prompt is purely
exploratory and makes no code or documentation change, do not bump the version.

## Commit And Push Expectations

When the user prompts for substantive implementation work, agents should prepare
clean commits and push code after the requested work is complete, unless the
user explicitly asks not to.

Before committing:

- inspect `git status --short`;
- verify only intended files are staged;
- run the most relevant tests/checks available in the environment;
- mention any checks that could not be run.

Do not tag releases without explicit approval.

## Roadmap And TODO

- Use `ROADMAP.md` for milestone intent and architecture.
- Use `TODO.md` for discrete implementation slices.
- Mark TODO items complete only after the work is actually implemented and
  verified.
- If new work is discovered, add it as a discrete checkbox under the appropriate
  slice rather than burying it in prose.
