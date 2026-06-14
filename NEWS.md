# floundeR 0.13.1

* Documented the Porkchop capabilities that remain out of scope for normal
  floundeR APIs, including read trimming, transformed FASTQ output, barcode
  demultiplexing, cDNA rescue, benchmarking, terminal dashboards, standalone
  HTML reports, custom ungoverned motif registries, and raw CLI passthrough.

# floundeR 0.13.0

* Added `library_preparation_report_card()` and
  `library_preparation_report_card_thresholds()` to convert Porkchop-derived
  kit, adapter/primer, barcode/flank, and cDNA primer evidence into the
  standard floundeR report-card schema for QC review and reporting.

# floundeR 0.12.0

* Added in-process Porkchop-backed R wrappers for library kit candidates,
  adapter/primer evidence, barcode/flank evidence, and cDNA primer-pair
  evidence. Until Porkchop is public, default source builds keep that engine
  unlinked and the wrappers fail with typed R conditions; local sibling-checkout
  builds can enable the `porkchop-integration` feature.
* Porkchop evidence preserves heuristic score terminology, kit provenance,
  lifecycle/support status, and validation limitations while keeping trimming
  and preprocessing out of the normal floundeR API.

# floundeR 0.11.4

* Added the Porkchop integration audit for curated library-preparation QC
  surfaces, preserving heuristic score terminology and excluding trimming-first
  APIs from the normal floundeR surface.

# floundeR 0.11.3

* Added an executable Synoptikon QC handoff vignette covering run QC evidence,
  ONT POD5 provenance, reserved BAM/library-preparation sections, and JSON
  payload writing.
* Updated the release-tarball check to build and check package vignettes.

# floundeR 0.11.2

* Added offline Synoptikon QC payload JSON-schema validation tests and fixed
  nested provenance metadata plus report-card check schema compatibility.

# floundeR 0.11.1

* Strengthened Synoptikon QC payload section handling so supplied flowcell,
  barcode, POD5, BAM, and library-preparation evidence contributes section
  status and is covered by offline tests.

# floundeR 0.11.0

* Added `as_synoptikon_qc()` and `write_synoptikon_qc()` to assemble existing
  floundeR QC evidence into the v1 Synoptikon payload contract and write it as
  JSON for downstream ingestion.

# floundeR 0.10.0

* Defined the first versioned Synoptikon QC payload contract and installed the
  v1 JSON schema for downstream ingestion, review, reporting, and audit
  handoff work.

# floundeR 0.9.2

* Documented Bamana capabilities that remain out of scope for floundeR so the
  BAM surface stays curated around nanopore QC, review, provenance, reporting,
  and synoptikon handoff rather than becoming a wholesale Bamana mirror.

# floundeR 0.9.1

* Documented the Bamana read-only versus transformation boundary, including
  which operations may feed QC evidence directly and which write-capable
  operations require explicit output paths and provenance manifests.

# floundeR 0.9.0

* Added `bam_qc_report_card()` and `bam_qc_report_card_thresholds()` to convert
  Bamana-derived summary, index, sort, tag, EOF, validation, and provenance
  evidence into the standard floundeR pass/warn/fail report-card schema.

# floundeR 0.8.1

* Documented and tested Bamana response-envelope semantics: payload-bearing
  QC failures remain structured evidence, while no-payload failures become
  typed R conditions carrying Bamana code, detail, and hint metadata.

# floundeR 0.8.0

* Added Rust-backed `bam_check_index()`, `bam_check_map()`,
  `bam_check_sort()`, and `bam_check_tag()` wrappers over Bamana's curated
  read-only report-card checks for index, mapping, sorting, and aux-tag
  evidence.

# floundeR 0.7.0

* Added Rust-backed `bam_verify()`, `bam_validate()`, and `bam_check_eof()`
  wrappers over Bamana's read-only verification, validation, and BGZF EOF
  evidence surfaces, with stable R-native return shapes and offline BAM
  fixture coverage.

# floundeR 0.6.0

* Added Rust-backed `bam_summary()` as the first curated Bamana integration,
  returning stable R-native QC tables for command status, evidence scope,
  header metadata, counts, fractions, MAPQ, mapping state, flags, references,
  index-derived evidence, and optional MAPQ histograms.

# floundeR 0.5.2

* Audited Bamana's public Rust command modules, JSON contracts, and governing
  docs, and recorded the curated read-only BAM QC/review surface that should be
  promoted into floundeR through in-process Rust bindings.

# floundeR 0.5.1

* Added an opt-in derived-POD5 demonstration workflow script that records ONT
  source metadata, read-only subdivision plans, and provenance manifests outside
  the repository without writing or committing derived POD5 files.

# floundeR 0.5.0

* Added read-only Rust-backed `pod5_subdivide_plan()` for deterministic POD5
  file-count, sample-label, elapsed-time, and read-count planning without
  writing derived POD5 files.

# floundeR 0.4.1

* Recorded upstream mocked-reader coverage for mixed-flow-cell POD5 folder
  aggregation and clarified that public R mixed-flow-cell tests require the
  future parser-backed `pod5-tools` reader.

# floundeR 0.4.0

* Added Rust-backed `pod5_folder_info()`, `pod5_manifest()`, and
  `pod5_compare()` wrappers for run-tree POD5 QC, versioned collection
  manifests, and operational handoff comparisons.

# floundeR 0.3.0

* Added Rust-backed `pod5_verify()` and `pod5_file_info()` wrappers that return
  R-native data frames while preserving the current `pod5-tools` fast
  extension/signature verification semantics and typed POD5 error conditions.

# floundeR 0.2.6

* Documented the no-commit policy for downloaded POD5 files and added
  repository ignore guardrails for local ONT open-data caches while preserving
  an explicit path for intentionally tiny test fixtures.

# floundeR 0.2.5

* Added no-network ONT Zymo fecal POD5 dataset and selected pass/fail example
  object helpers so documentation, downloads, and future reports share stable
  source metadata.

# floundeR 0.2.4

* Added `ont_open_data_fetch()` for explicit single-object ONT open-data
  downloads into caller-controlled cache directories with returned provenance
  metadata.

# floundeR 0.2.3

* Added `ont_open_data_list()` for explicit anonymous metadata listings of the
  selected ONT Zymo fecal POD5 S3 prefix, including an opt-in documentation
  example that does not download POD5 files.

# floundeR 0.2.2

* Added an executable offline `pod5_find()` example that demonstrates
  in-process Rust-backed folder discovery without invoking an external CLI or
  requiring network/downloaded POD5 data.

# floundeR 0.2.1

* Documented the current `R CMD check` compiled-code warning introduced by
  linking the Rust `std`/`pod5-tools` stack, including the release posture for
  GitHub builds and the requirement to resolve or obtain policy acceptance
  before Bioconductor/CRAN-facing release work.

# floundeR 0.2.0

* Added `pod5_find()` as the first curated POD5 Rust-backed R API. The function
  calls the `../pod5-tools` discovery library in-process and returns an R data
  frame with POD5-containing folders, file counts, byte totals, and modification
  windows.

# floundeR 0.1.10

* Raised the Rust toolchain floor to 1.85 and moved the embedded crate to Rust
  edition 2024 so floundeR can depend on the current `../pod5-tools` library.

# floundeR 0.1.9

* Recorded the POD5 discovery preflight against `../pod5-tools`, confirming
  that the existing `find_pod5_directories()` library API is the intended
  curated source for `pod5_find()`.

# floundeR 0.1.8

* Documented source-install requirements for macOS, Linux, Windows, Docker, and
  CI builds of the embedded Rust extension, and made CI set up Rust explicitly.

# floundeR 0.1.7

* Added public R wrappers for compiled Rust capability checks, typed
  `floundeR_rust_unavailable` conditions, and a downstream test skip helper.

# floundeR 0.1.6

* Added an internal Rust capability function callable from R, proving the
  `extendr` scaffold can return structured R data through registered native
  symbols.

# floundeR 0.1.5

* Added the initial `extendr` scaffold under `src/rust` with R package build
  hooks and development-container Rust tooling.

# floundeR 0.1.4

* Recorded the Slice 6 cross-repository Rust preflight, confirming that the
  first local `extendr` scaffold does not require adjacent changes to
  `pod5-tools`, `bamana`, or `porkchop`.

# floundeR 0.1.3

* Recorded the Rust-in-R architecture decision for an embedded `extendr` crate
  with curated path dependencies on the POD5, BAM/BGZF/FASTQ, Porkchop, and
  optional private Grammateus engine libraries.

# floundeR 0.1.2

* Added `qc_report_card()` with schema-versioned pass/warn/fail checks and
  documented default thresholds for sequencing-summary-derived QC evidence.

# floundeR 0.1.1

* Added tidy, schema-versioned QC contracts for yield over time, read length
  distribution, quality distribution, channel density, and barcode composition.

# floundeR 0.1.0

* Added `qc_run_summary()` as the first versioned core QC contract for
  nanopore run-level summaries from normalised sequencing-summary data.

# floundeR 0.0.24

* Reconciled the current import surface by declaring narrow namespace imports
  for active FASTA, FASTQ, GenBank, plotting, and ONT open-data dependencies.
* Removed stale RMarkdown-era packages from the development dependency
  bootstrap list now that legacy vignettes are excluded from package checks.
* Documented the remaining legacy `%>%` and `BamFile` exports so release-style
  package checks can complete without warnings.

# floundeR 0.0.23

* Added a release-style package check helper that builds and checks a source
  tarball outside the repository, avoiding diagnostics from development-only
  files and local check output.
* Corrected source-build ignore rules so repository governance files,
  development scripts, `.github`, and `.dockerignore` are excluded from package
  tarballs.

# floundeR 0.0.22

* Moved legacy RMarkdown vignettes to repo-only historical reference material
  while executable report workflows move to the Grammateus migration.

# floundeR 0.0.21

* Corrected POD5 manifest fixture checksum values so they remain character
  SHA-256-like hex strings under package checks and schema tests.

# floundeR 0.0.20

* Added a small installed BLAST fixture and restored `Blast$count()` behavior
  so the BLAST example and test no longer depend on a missing large UniRef
  output file.

# floundeR 0.0.19

* Added conservative sequencing-summary column alias handling for current
  Dorado and MinKNOW/POD5-era summary shapes.
* Allowed reduced Dorado summary files without `passes_filtering` to import as
  partial QC inputs with an explicit warning and `NA` logical values.

# floundeR 0.0.18

* Replaced generated `eval(parse())` sequencing-summary column specifications
  with structured `readr::cols_only()` collector construction.
* Added regression coverage for legacy Guppy-like and Dorado-like sequencing
  summary fixtures.

# floundeR 0.0.17

* Added a Docker-backed development build path and a dependency audit script for
  assessing declared CRAN/Bioconductor package versions and support channels in
  a clean container.
* Recorded the current R 4.6.0/Bioconductor 3.23 dependency snapshot and the
  package-quality failures exposed once dependencies are available.

# floundeR 0.0.16

* Added fixture-backed schema tests for the initial sequencing-summary,
  barcode-summary, and POD5-manifest QC table contracts.

# floundeR 0.0.15

* Added a test helper and development guidance so optional compiled Rust-backed
  tests skip cleanly until in-process bindings are available.

# floundeR 0.0.14

* Added a test helper and development guidance that keep network-dependent ONT
  open-data/S3 tests opt-in through `FLOUNDER_RUN_NETWORK_TESTS=true` and
  skipped by default on CRAN-like checks.

# floundeR 0.0.13

* Added small offline test fixtures for sequencing summaries, barcoding
  summaries, FASTQ, FASTA, and POD5 metadata manifests without committing
  binary POD5 files.

# floundeR 0.0.12

* Added a dry-run-first R dependency bootstrap script and development notes for
  local checks after FAST5 retirement.
* Configured the package for `testthat` edition 3.

# floundeR 0.0.11

* Retired the active FAST5 API surface by removing the `Fast5` export,
  implementation, fixtures, tests, vignette, generated reference page, and
  `rhdf5` dependency. Raw-signal QC now points toward POD5-era integration.

# floundeR 0.0.10

* Added a repository-level governance check that enforces the open-source
  boundary for floundeR, pod5-tools, bamana, and porkchop while keeping
  Grammateus private, optional, and distributed through prebuilt runtime assets.

# floundeR 0.0.9

* Added Porkchop to the curated open-source Rust engine plan for adapter,
  primer, barcode, kit-registry, and cDNA/library-preparation QC evidence.

# floundeR 0.0.8

* Documented the distribution boundary: `floundeR`, `pod5-tools`, and `bamana`
  stay open-source, while Grammateus remains private and is consumed through
  optional prebuilt runtime/reporting artifacts.

# floundeR 0.0.7

* Added an explicit Grammateus reporting interface plan for dropping R plots
  into governed HTML/PDF reports and for controlled R plot generation from
  Grammateus report specifications.

# floundeR 0.0.6

* Started the reboot planning baseline for contemporary nanopore QC and review.
* Documented the FAST5 retirement direction, selective Rust-in-R integration
  with `pod5-tools`, `bamana`, and `grammateus`, ONT POD5 example-data policy,
  and synoptikon reporting ambitions.
* Reconciled package metadata and declared currently referenced dependencies so
  package checks can progress beyond initial DESCRIPTION validation.

# floundeR 0.0.5

# floundeR 0.0.4

* inclusion of a basic Rsamtools derived FASTA R6 class (lodestar dependency)
* starting to ensure that methods are also available through documented
  magrittr %>% pipes for legibility of code
* robust Angenieux decorations and demonstration of more comprehensive plots
* integrating Angenieux with ggsave for cleaner Rmarkdown usage

# floundeR 0.0.3

* pkgdown documentation for the github pages
* creation of angenieux prototype; code linking and POC
* creation of the R6 object called FlowCell
* Sequencing summary push ...
* simplified FAST5.is_multi_F5 logic; there were some issues with different
  version of the hdf5 api - tested on Windows and macOS. 

# floundeR 0.0.2

* starting to populate content
* inclusion of FAST5 parsing content, vignette and test data

# floundeR 0.0.1

* tagged a first version

# floundeR 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
