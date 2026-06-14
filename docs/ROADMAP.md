# floundeR Revival Roadmap

This roadmap reboots `floundeR` as a contemporary R package for Oxford
Nanopore QC, with POD5-aware raw-data inspection backed by Rust code
callable inside R. The package should not shell out to a Rust CLI for
normal operation. Rust should be used as an in-process implementation
layer exposed through R functions and R6 classes. The reboot is not
exclusively focused on POD5: basecalled and aligned outputs, especially
BAM, remain first-class QC inputs. Library-preparation evidence such as
adapters, primers, barcodes, and kit identity also belongs in the
QC/review surface.

## Code Review Findings

### P0: Package Check Fails Before Code Is Loaded

`R CMD check --no-manual --no-build-vignettes .` currently fails on this
machine with:

``` text
Required fields missing or empty:
  'Author' 'Maintainer'
```

The `DESCRIPTION` uses `Authors@R`, which should normally be sufficient,
but the R 4.6 check environment is still requiring generated legacy
fields. The revival needs a first packaging slice that makes
`R CMD check` pass metadata validation in the target toolchain.

### P0: FAST5 Is Still a First-Class API Surface

FAST5 remains exported as `Fast5`, documented, tested, included in
fixtures, and backed by `rhdf5`. This conflicts with the project
direction. The removal is not just deleting `R/fast5.R`; it also
requires coordinated changes in `NAMESPACE`, `DESCRIPTION`, `man/`,
`vignettes/`, `docs/`, fixtures, NEWS, and tests.

### P1: FAST5 Implementation Is Not Stable Enough to Preserve

The multi-FAST5 path samples a random top-level read key to infer file
shape. That means metadata returned by the class can vary by run and
test coverage may pass for the wrong reason when files contain
heterogeneous metadata. Since the feature is being retired, do not
invest in repairing this logic except for a temporary deprecation shim
if needed.

### P1: Sequencing Summary Parsing Uses Generated R Code

`SequencingSummary` builds a `readr::cols_only(...)` expression as text
and then calls `eval(parse(...))`. This should be replaced with a
structured column-spec builder. The current approach makes column
parsing harder to test, harder to lint, and brittle when Dorado or
MinKNOW column names evolve.

### P1: Declared Dependencies Need Reconciliation

The code and generated namespace reference packages that are not
declared in `DESCRIPTION`, including `GenomeInfoDb`, `GenomicRanges`,
`IRanges`, `BiocGenerics`, `RColorBrewer`, `magick`, `stringdist`, and
`tidyselect`. Conversely, `rhdf5` should be removed with FAST5.
Dependency cleanup must be treated as part of package revival, not
cosmetic maintenance.

### P2: Documentation Is Generated But Stale

The repository carries generated `man/` and `docs/` artifacts that still
advertise FAST5. Until regenerated, the public story will be internally
contradictory. The reboot should define a generated-doc policy: either
commit pkgdown output intentionally after release-oriented changes, or
remove generated site files from normal development.

### P2: Synoptikon/Mnemosyne Integration Needs Explicit Contracts

The package has useful QC concepts, but it does not yet expose stable
machine-readable result contracts for `../mnemosyne`/synoptikon. The
reboot should define return schemas for run summaries, flowcell density,
barcode yield, POD5 folder integrity, manifest comparison, and report
cards.

## Positioning

`floundeR` should become the R-native QC and review layer for nanopore
sequencing runs:

- fast import of sequencing summaries, BAM/FASTQ summaries, and POD5
  metadata;
- run-level QC objects that return tidy data frames and
  publication-ready plots;
- selective in-process Rust acceleration for POD5 discovery,
  verification, manifesting, subsetting support, and provenance through
  `../pod5-tools`;
- selective in-process Rust acceleration for BAM/BGZF/FASTQ QC evidence
  through `../bamana`;
- selective in-process Rust acceleration for adapter, primer, barcode,
  kit, and cDNA/library-preparation evidence through `../porkchop`;
- report rendering through Grammateus semantic report contracts rather
  than RMarkdown-first documents;
- stable JSON/data-frame contracts for synoptikon ingestion;
- a package API that feels native to R users while handling
  operationally large long-read datasets.

The package should avoid duplicating official ONT tools where those
tools are already excellent. Its differentiator should be operational QC
and integration: run trust, provenance, metadata consistency,
storage-tree audits, flowcell health, barcode performance, and compact
handoff reports.

This is not a plan to shoehorn every feature from `../pod5-tools` or
`../bamana` into an R package. Those projects should remain strong
specialist toolkits. `floundeR` should use them wisely as Rust
implementation engines for the subset of behavior that directly improves
nanopore sequence-data QC, review, reporting, and synoptikon handoff.

Distribution boundary: `floundeR`, `pod5-tools`, `bamana`, and
`porkchop` are open-source. Grammateus remains private. `floundeR` must
therefore support optional prebuilt Grammateus runtime/reporting
artifacts for HTML/PDF generation without making private Grammateus
source mandatory for public installation, checking, or core QC use.

## Rust-In-R Integration Strategy

The right target is an in-process Rust layer, not a CLI wrapper.

Recommended shape:

- add a Rust extension crate under `src/rust/` or equivalent R package
  conventions;
- use `extendr` to expose Rust functions to R;
- depend on reusable logic from `../pod5-tools` as a Rust library
  dependency, either by path during local development or by a released
  git/crates version;
- depend on reusable logic from `../bamana` as a Rust library dependency
  for BAM, BGZF, FASTQ, ingest, sampling, summary, validation, and
  forensic QC surfaces;
- depend on reusable logic from `../porkchop` as a Rust library
  dependency for adapter, primer, barcode, kit-registry, cDNA, and
  library-preparation QC evidence;
- keep CLI-only concerns in `pod5-tools` out of the R binding surface;
- expose R functions only where they serve QC/review workflows, such as
  [`pod5_find()`](https://sagrudd.github.io/floundeR/reference/pod5_find.md),
  [`pod5_file_info()`](https://sagrudd.github.io/floundeR/reference/pod5_file_info.md),
  [`pod5_verify()`](https://sagrudd.github.io/floundeR/reference/pod5_verify.md),
  [`pod5_folder_info()`](https://sagrudd.github.io/floundeR/reference/pod5_folder_info.md),
  selected manifest/provenance helpers, and demonstration subsetting
  helpers;
- expose R functions only where they serve QC/review workflows, such as
  [`bam_summary()`](https://sagrudd.github.io/floundeR/reference/bam_summary.md),
  [`bam_verify()`](https://sagrudd.github.io/floundeR/reference/bam_verify.md),
  [`bam_validate()`](https://sagrudd.github.io/floundeR/reference/bam_validate.md),
  selected EOF/index/map/sort and tag checks, and BAM-aware QC
  report-card builders;
- expose R functions only where they serve QC/review workflows, such as
  `library_screen()`, `kit_candidates()`, `adapter_primer_evidence()`,
  `barcode_evidence()`, and `cdna_primer_evidence()`;
- return tibbles/lists with stable column names, not rendered TSV/JSON
  strings;
- add R6 wrappers only where they improve workflow composition.

The Rust library currently contains many useful contracts, but it still
includes CLI parsing and rendered-output paths. The reboot should
extract or harden the pure library APIs so R can call typed functions
directly.

## Bamana BAM Processing Strategy

`floundeR` should link BAM processing functionality from `../bamana`
rather than reimplement BAM parsing in R. Bamana is a Rust toolkit for
verification, quality control, inspection, and transformation of BAM and
adjacent formats. Its local charter states that performance-critical
BAM, BGZF, FASTQ, sampling, ingest, and forensic hot paths must remain
Bamana-native, with external crates kept to compatibility/oracle roles.

Target shape:

- integrate Bamana through in-process Rust functions callable from R;
- do not shell out to the `bamana` CLI for normal package APIs;
- preserve Bamana’s JSON/schema contract semantics while converting
  outputs to R-native data frames/lists and conditions;
- use Bamana for alignment-level QC: mapping summaries, MAPQ, flags,
  duplicate and QC-fail fractions, index state, sort state, EOF
  evidence, tag presence, validation findings, checksums, subsampling,
  filtering, and forensic provenance signals;
- avoid exposing Bamana’s entire transformation surface unless a
  function has a clear QC/review purpose in `floundeR`;
- keep header-only mutation, record-level read-group annotation,
  biological duplicate marking, collection-duplication remediation, and
  provenance inspection conceptually distinct, matching the Bamana
  charter;
- allow `../bamana` modifications where needed to expose clean library
  APIs for `floundeR`, while preserving Bamana’s charter and
  documentation contracts.

No bamana-specific canonical files were found under `../mnemosyne-docs`
during the initial review. Until such files exist, treat
`../bamana/CHARTER.md`, `../bamana/docs/project-charter.md`, and
Bamana’s public JSON/schema documentation as the governing product
boundary.

## Porkchop Library-Preparation Strategy

`floundeR` should link selected functionality from `../porkchop` for
library-preparation QC. Porkchop is an open-source Rust toolkit for ONT
kit inspection, adapter/barcode screening, trimming, cDNA workflow
support, and benchmarking. The floundeR integration should focus on
evidence and review, not on turning floundeR into a general
trimming/preprocessing frontend.

Target shape:

- integrate Porkchop through in-process Rust functions callable from R;
- do not shell out to the `porkchop` CLI for normal package APIs;
- use Porkchop for kit-registry metadata, adapter/primer/barcode motif
  evidence, cDNA primer-pair evidence, and library chemistry mismatch
  signals;
- preserve Porkchop’s explicit language that screening scores are
  heuristic evidence scores, not calibrated probabilities;
- preserve kit provenance, lifecycle status, support level, and
  validation limitations in floundeR report outputs;
- allow `../porkchop` modifications where needed to expose clean library
  APIs for `floundeR`, while preserving Porkchop’s AGENTS, roadmap,
  Sphinx docs, JSON schema contracts, and kit provenance rules.

No porkchop-specific canonical files were found under
`../mnemosyne-docs` during the initial review. Until such files exist,
treat `../porkchop/AGENTS.md`, `../porkchop/ROADMAP.md`,
`../porkchop/README.md`, and Porkchop’s Sphinx output
contract/provenance documentation as the governing product boundary.

## Grammateus Reporting Strategy

Report rendering should move away from RMarkdown as the primary report
mechanism. RMarkdown can remain as a temporary compatibility path while
the package adopts Grammateus fully.

The retained RMarkdown boundary is now explicit: RMarkdown is allowed
for package vignettes, migrated tutorials, and transitional examples,
but not for production QC report generation, Mnemosyne Biosciences
branded technical reports, or Synoptikon trusted-report lifecycle
artifacts. The operational policy is documented in
`LEGACY_REPORTING.md`.

Target shape:

- use Grammateus as an in-process Rust report-rendering and
  report-contract layer inside R;
- do not shell out to Grammateus command-line binaries for normal report
  rendering;
- support both directions of control: R can pass existing `ggplot2` plot
  artifacts into Grammateus reports, and Grammateus can call controlled
  R plot generation from semantic report specifications;
- represent QC tables, figures, plots, images, report cards, provenance,
  and report lifecycle metadata as Grammateus semantic elements;
- keep Mnemosyne Biosciences branding, typography, and technical-report
  styling in shared Grammateus templates/themes rather than bespoke
  RMarkdown CSS;
- produce HTML/PDF/report-manifest artifacts through Grammateus
  contracts;
- carry provenance for source data, software version, production
  timestamp, plot/table hashes, and report payload hashes;
- align trusted-report output with Synoptikon/Mnemosyne lifecycle and
  governance-signing contracts.
- document the R-facing interface in `REPORTING_INTERFACE.md`.

The local `../grammateus` checkout is source-complete and exports report
element, rendering, provenance, trusted-report, governance, template,
and registry APIs from its Rust library. Canonical product documentation
lives in `../mnemosyne-docs/products/grammateus`. A reporting
integration slice should depend on the Grammateus library directly and
add any missing public library APIs deliberately, with tests on both
sides.

For open-source distribution, that development-time private checkout
must be converted into a prebuilt optional runtime interface. The public
floundeR package should detect, validate, and use a Grammateus runtime
bundle when available, but it must install and provide core QC
functionality without that bundle.

The distribution model is documented in `DISTRIBUTION.md`.

Cross-repository changes are in scope for this reboot. `../pod5-tools`,
`../bamana`, `../porkchop`, and `../grammateus` may be modified where
needed to provide clean library APIs for `floundeR`, but those changes
must preserve the relevant product charters, requirements, architecture
decisions, and canonical contracts maintained in `../mnemosyne-docs` or
in the adjacent repository when no canonical docs exist there yet.

## External Data Strategy

The original package used S3-hosted example data. The reboot should keep
that strength, but update it for POD5-era workflows and avoid committing
large raw files.

Canonical POD5 integration dataset:

``` text
s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/
```

This public ONT Zymo fecal metagenomics run was verified on 2026-06-13
through anonymous S3 listing in `eu-west-1`. The prefix contains 9107
POD5 objects: 8290 pass files and 817 fail files, totalling about 2.68
TB. Individual files are large, so package tests must use small local
fixtures and keep ONT downloads opt-in.

Fixed example-file decision:

- primary pass example:
  `zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_pass_279c9095_68316534_8289.pod5`
  at `47077200` bytes;
- fail-state example:
  `zymo_fecal_2025.05/raw/PAU85136/pod5/PAU85136_fail_279c9095_68316534_0.pod5`
  at `163007608` bytes.

The pass example should be the default real-data demonstration. The fail
example should be downloaded only for examples that need failure-state
QC. `../pod5-tools` can then derive smaller POD5 subsets or splits from
these selected source files for demonstrations, tests, reports, and
synoptikon demos. Those derived artifacts must retain provenance linking
back to the ONT S3 source object and the `pod5-tools` parameters/version
used to create them.

Target shape:

- use `aws.s3` for anonymous bucket listing and controlled downloads;
- wrap ONT open-data access in explicit helper functions rather than
  embedding S3 logic throughout vignettes;
- cache downloaded examples outside the repository;
- use the selected Zymo fecal POD5 files for integration tests, report
  examples, and performance benchmarks;
- use `pod5-tools` to create smaller derived POD5 subsets/splits when
  full source objects are unnecessarily large for a demonstration;
- record source object keys, sizes, timestamps, and checksums/manifests
  in QC provenance.

## Milestone 0: Stabilise The Package Baseline

Goal: make `floundeR` installable, checkable, and honest about current
support.

Deliverables:

- repair package metadata so `R CMD check` reaches code checks;
- document the reboot direction in README and NEWS;
- add CI for R package checks on the intended R versions;
- decide whether `docs/` remains committed or is release-generated only;
- capture the current exported API and deprecation plan.
- record cross-repository constraints for `../pod5-tools`,
  `../grammateus`, and `../mnemosyne-docs`.
- record cross-repository constraints for `../porkchop`.
- record the open-source/private distribution boundary for floundeR,
  pod5-tools, bamana, and Grammateus.

Acceptance criteria:

- `R CMD check --no-manual --no-build-vignettes .` gets past DESCRIPTION
  validation;
- dependency declarations match the R code;
- CI has at least one R package check job;
- FAST5 retirement is visible to users before deletion lands.

## Milestone 1: Retire FAST5

Goal: remove FAST5 capabilities and the `rhdf5` dependency without
leaving broken documentation, examples, or tests.

Deliverables:

- remove `R/fast5.R`, `Fast5` export, FAST5 fixtures, FAST5 tests, and
  FAST5 vignettes;
- remove `rhdf5` and FAST5 installation instructions;
- clean generated documentation and NEWS;
- optionally add a short deprecation note explaining that POD5 is the
  supported raw-signal direction.

Acceptance criteria:

- no package code imports `rhdf5`;
- `rg -i "fast5|rhdf5"` only finds historical NEWS or explicit migration
  notes;
- tests do not require FAST5 fixtures;
- package check remains green.

## Milestone 2: Modernise Core QC APIs

Goal: keep and sharpen the R-native parts that make the package
valuable.

Deliverables:

- replace dynamic parsing in sequencing summary import with structured
  column specifications;
- support Dorado/MinKNOW-era sequencing summary column variants;
- define a single run-level QC object or function family for summary
  metrics;
- return tidy data frames for yield, quality, temporal progression,
  channel density, and barcode composition;
- add snapshot or schema tests for all user-facing returned data.

Acceptance criteria:

- sequencing summary import is deterministic and does not use
  `eval(parse())`;
- missing optional columns produce clear warnings and partial results;
- public QC returns have documented columns;
- core tests run without network access.

## Milestone 3: Add In-Process POD5 Rust Binding

Goal: expose POD5 operational metadata to R through compiled Rust code.

Deliverables:

- add an `extendr` build scaffold;
- expose typed Rust functions for POD5 discovery and validation;
- connect to `pod5-tools` library APIs by path during development;
- convert Rust structs to R data frames/lists;
- add installation guidance for users building from source;
- add tests that skip cleanly when Rust tooling is unavailable.
- add opt-in integration examples against the ONT Zymo fecal POD5
  dataset.

Acceptance criteria:

- R users call `pod5_find(path)` and receive a data frame;
- R users call `pod5_verify(path)` and receive a structured result;
- no normal R API shells out to `pod5-tools`;
- error classes are predictable and test-covered.

## Milestone 4: Add In-Process BAM/Bamana Rust Binding

Goal: expose BAM operational QC to R through compiled Rust code from
`../bamana`.

Deliverables:

- identify Bamana library functions and CLI-only paths that need
  promotion to reusable APIs for QC/review purposes;
- expose typed Rust functions for BAM summary, verification, validation,
  EOF, index, map, sort, tag, checksum, and bounded/full scan QC;
- convert Bamana response envelopes into R data frames/lists while
  preserving command semantics and warnings;
- define BAM QC report-card checks for mapping fraction, unmapped
  fraction, duplicate fraction, QC-fail fraction, MAPQ zero burden,
  missing tags, stale indices, sort/index mismatch, EOF absence,
  validation findings, and provenance/forensic anomalies;
- add tests using small BAM fixtures and skip rules for optional large
  integration examples.

Acceptance criteria:

- R users call `bam_summary(path)` and receive a stable structured
  result;
- R users call `bam_validate(path)` or `bam_verify(path)` and receive
  actionable findings;
- no normal R API shells out to `bamana`;
- Bamana charter distinctions are preserved in R names and docs.
- the R API is curated around nanopore QC and review rather than a
  wholesale mirror of Bamana.

## Milestone 5: Synoptikon QC Contracts

Goal: make `floundeR` useful as the QC ingestion layer for synoptikon
via `../mnemosyne`.

Deliverables:

- define
  [`qc_run_summary()`](https://sagrudd.github.io/floundeR/reference/qc_run_summary.md)
  returning a compact, versioned list/data frame;
- define
  [`qc_report_card()`](https://sagrudd.github.io/floundeR/reference/qc_report_card.md)
  for pass/warn/fail operational checks;
- define schemas for flowcell, temporal, barcode, POD5 folder, BAM
  alignment, BAM validation, BAM index/sort/map evidence,
  library-preparation evidence, and manifest metrics;
- add JSON export helpers for synoptikon ingestion;
- add fixtures that represent realistic pass, warning, and failure
  scenarios.
- include provenance fields for external ONT open-data object keys when
  reports are generated from downloaded example data.

Acceptance criteria:

- synoptikon can ingest one run QC payload without bespoke parsing;
- schema versions are included in exported QC payloads;
- validation tests catch accidental column or type drift.

## Milestone 6: Grammateus Reporting And Visual Identity

Goal: make the package something people want to show, not just install.

Deliverables:

- refresh README around current nanopore workflows;
- add a modern pkgdown site with examples focused on QC questions;
- add a run QC vignette built around real sequencing-summary/POD5
  metadata;
- add Grammateus-backed technical report rendering for run QC;
- define optional Grammateus runtime discovery, installation,
  validation, and compatibility checks for GitHub-distributed prebuilt
  artifacts;
- add a documented R/Grammateus plot interface covering R-to-Grammateus
  artifact handoff and Grammateus-to-R controlled plot generation;
- define Grammateus semantic elements for run metadata tables, yield
  plots, flowcell density figures, barcode summaries, POD5 integrity
  tables, report card findings, BAM alignment summaries, BAM
  validation/index/sort evidence, library-preparation evidence,
  kit-candidate summaries, adapter/primer/barcode findings, methods,
  limitations, and appendices;
- encode Mnemosyne Biosciences branding and styling through Grammateus
  templates/themes;
- reuse existing Grammateus public exports such as `ReportTable`,
  `ReportFigure`, `ReportPlot`, `ReportElementProvenance`,
  trusted-report manifest types, and HTML/PDF rendering functions
  wherever possible;
- emit report manifests and provenance metadata suitable for
  trusted-report lifecycle handling;
- refine plots for flowcell density, cumulative yield, quality
  distribution, barcode balance, acquisition gaps, and POD5 integrity;
- add a one-command HTML/PDF report path backed by Grammateus.

Acceptance criteria:

- the first README example performs a useful QC workflow;
- vignettes are current and executable or intentionally precomputed;
- plot functions are tested for return classes and minimal data
  requirements;
- docs no longer mention retired FAST5 workflows except in migration
  notes.
- rendered reports are generated through in-process Rust/Grammateus
  APIs, not RMarkdown-only documents or external CLI calls;
- reports include stable element identifiers, captions, alt text where
  applicable, provenance, and Mnemosyne Biosciences branding.
- existing R plots can be dropped into reports through a documented
  helper, and semantic Grammateus plot specs can trigger controlled
  R/ggplot2 rendering.
- core floundeR QC functionality remains usable without private
  Grammateus source or runtime assets.

## Milestone 6: Release Discipline

Goal: move from revival to maintained package.

Deliverables:

- semantic versioning policy in `AGENTS.md` and README;
- `NEWS.md` updated for every substantive change;
- release checklist covering checks, docs, examples, and synoptikon
  contract compatibility;
- branch/commit/PR habits suitable for automation-assisted development.

Acceptance criteria:

- each substantive agent prompt bumps patch or minor version as
  appropriate;
- package checks and tests are run before handoff;
- releases have coherent NEWS entries and git tags after explicit
  approval.

## What To Add

High-value additions for a package people should talk about:

- [`qc_run_summary()`](https://sagrudd.github.io/floundeR/reference/qc_run_summary.md):
  one call that produces the core run metrics needed by humans and
  synoptikon.
- `pod5_*()` functions backed by in-process Rust: discovery, verify,
  file info, folder info, manifest, compare, and subdivision planning.
- Dorado-era compatibility: ingest current sequencing-summary and
  barcode summary outputs without assuming old Guppy names only.
- Report cards: pass/warn/fail checks for low yield, poor quality,
  channel dropouts, barcode imbalance, mixed flow cells, broken POD5
  files, acquisition gaps, and missing metadata.
- Stable schemas: version every exported QC payload.
- Modern docs: a README that answers “is this run healthy?” within the
  first minute.
- Small realistic fixtures: enough data to test behavior without
  shipping large raw datasets.
- Synoptikon exporter: write exactly the QC payload expected by
  `../mnemosyne`.
