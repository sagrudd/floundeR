# floundeR Revival TODO

Work from top to bottom unless a later slice is explicitly pulled forward.
Each slice should be small enough to commit independently. Substantive work
should include a semantic version bump and `NEWS.md` entry.

## Slice 1: Package Baseline

- [x] Fix `DESCRIPTION` metadata so `R CMD check` passes initial validation on
      the target R toolchain.
- [x] Reconcile `DESCRIPTION` dependencies with packages used in `R/`.
- [x] Confirm `aws.s3` remains an intentional dependency for ONT open-data
      access, or replace it with a documented alternative before removing it.
- [x] Identify the relevant `../mnemosyne-docs` charters/contracts for
      `floundeR`, `pod5-tools`, `bamana`, `grammateus`, synoptikon, and
      Mnemosyne QC reporting before cross-repository implementation begins.
- [x] If no canonical `../mnemosyne-docs` entry exists for `bamana`, record
      `../bamana/CHARTER.md` and `../bamana/docs/project-charter.md` as the
      governing boundary for BAM integration until canonical docs are added.
- [x] If no canonical `../mnemosyne-docs` entry exists for `porkchop`, record
      `../porkchop/AGENTS.md`, `../porkchop/ROADMAP.md`, `../porkchop/README.md`,
      and relevant Sphinx docs as the governing boundary for
      library-preparation QC integration until canonical docs are added.
- [x] Decide and document whether generated `docs/` are committed during normal
      development or only during releases.
- [x] Add or update CI for R package check.
- [x] Record the FAST5 retirement direction in `NEWS.md`.
- [x] Record and enforce the distribution boundary: floundeR, pod5-tools,
      bamana, and porkchop are open-source; Grammateus remains private and
      optional from the public package's perspective.
- [x] Run `R CMD check --no-manual --no-build-vignettes .` and capture any
      remaining failures as follow-up tasks.

Baseline check result on 2026-06-13: `R CMD check` now passes DESCRIPTION and
namespace validation, then stops because the local R library lacks required
imports/suggested packages. The FAST5 retirement slice removed `rhdf5`; the
remaining missing packages are CRAN and Bioconductor dependencies needed by the
current non-FAST5 package surface.

## Slice 2: FAST5 Retirement

- [x] Remove `Fast5` from `NAMESPACE`.
- [x] Remove `R/fast5.R`.
- [x] Remove `rhdf5` from `DESCRIPTION`.
- [x] Remove FAST5 fixtures from `inst/extdata/`.
- [x] Remove `tests/testthat/test-fast5.R`.
- [x] Remove or replace `vignettes/03_fast5_files.Rmd`.
- [x] Remove generated `man/Fast5.Rd`.
- [x] Regenerate or clean pkgdown/docs references to FAST5.
- [x] Update README installation instructions to remove VBZ/rhdf5 guidance.
- [x] Verify `rg -i "fast5|rhdf5"` only finds accepted migration/history notes.

## Slice 3: Test Harness Revival

- [x] Add a dependency bootstrap note or script for local development,
      including CRAN and Bioconductor dependencies needed after FAST5
      retirement is complete.
- [x] Add a Docker-backed development build and dependency audit path for the
      current R package dependency surface.
- [x] Record current container dependency versions and support posture.
- [x] Ensure `testthat` is installed/declared with a current edition.
- [x] Add helper fixtures for sequencing summary, barcode summary, FASTQ, FASTA,
      and minimal POD5 metadata tests.
- [x] Ensure network-dependent ONT open-data tests are opt-in and skipped by
      default on CRAN-like/local package checks.
- [x] Make tests skip cleanly when optional compiled Rust support is absent.
- [x] Add schema tests for public QC return values.
- [x] Add a local command note for running focused tests and full checks.

## Slice 4: Sequencing Summary Modernisation

- [x] Replace `eval(parse())` column specification construction with structured
      `readr::cols_only()` construction.
- [x] Support current Dorado/MinKNOW sequencing summary naming variants.
- [x] Return clear partial-result warnings for absent optional columns.
- [x] Add tests for old Guppy-like and current Dorado-like summaries.
- [x] Correct documentation that describes sequencing-summary inputs as FAST5.

## Slice 4A: Check Debt Exposed By Container Build

- [x] Repair `Blast` example/test fixture resolution under `R CMD check`.
- [x] Correct POD5 manifest fixture checksum lengths.
- [x] Decide whether legacy vignettes should execute during development checks
      or be retired behind the Grammateus reporting migration.
- [x] Run release-style checks from an `R CMD build` tarball to avoid
      source-tree hidden-file noise.
- [x] Reconcile unused `Imports` after the legacy API surface is narrowed.

## Slice 5: Core QC Contracts

- [x] Define `qc_run_summary()` and its schema version.
- [x] Define tidy outputs for yield over time, read length distribution,
      quality distribution, channel density, and barcode composition.
- [x] Add `qc_report_card()` with pass/warn/fail checks.
- [x] Add tests for report-card thresholds and missing data.
- [x] Document return columns and threshold defaults.

## Slice 6: In-Process Rust Scaffold

- [x] Choose final Rust-in-R structure: embedded crate in `floundeR`, path
      dependencies on `../pod5-tools`, `../bamana`, `../porkchop`, and
      `../grammateus`, or split support packages.
- [x] Define a curation rule for Rust-backed APIs: include only functionality
      that directly supports nanopore QC, review, reporting, provenance, or
      synoptikon handoff.
- [x] If `../pod5-tools` changes are needed, verify they preserve the relevant
      `../mnemosyne-docs` charter and contracts before implementation.
- [x] If `../bamana` changes are needed, verify they preserve the Bamana
      charter, public JSON/schema contracts, and relevant docs before
      implementation.
- [x] If `../porkchop` changes are needed, verify they preserve Porkchop's
      AGENTS, roadmap, Sphinx docs, JSON schema contracts, kit provenance
      rules, and score-terminology requirements before implementation.
- [x] Add an `extendr` scaffold and build configuration.
- [x] Add a minimal Rust function callable from R.
- [x] Add R wrappers with typed errors and test skips when Rust tooling is
      unavailable.
- [x] Document source-install requirements for macOS, Linux, and CI.

## Slice 7: POD5 Discovery Binding

- [x] Promote only the `pod5-tools` logic needed for floundeR QC/review into
      pure library APIs if needed.
- [x] Align floundeR's Rust toolchain floor with `../pod5-tools` before adding
      the source-package-safe Rust dependency for `pod5_find()`.
- [x] Expose `pod5_find(path)` as an R function returning a data frame.
- [x] Include path, file count, byte total, oldest mtime, and newest mtime.
- [x] Add tests for empty folders, nested folders, and mixed file extensions.
- [ ] Resolve or formally document the release-check compiled-code warning from
      linking the Rust `std`/`pod5-tools` stack before Bioconductor/CRAN-facing
      release work.
- [ ] Add examples that do not call any external CLI process.
- [ ] Add an opt-in example that lists the ONT Zymo fecal POD5 S3 prefix after
      explicit user consent/network availability checks.

## Slice 7A: ONT Open-Data Helpers

- [ ] Implement `ont_open_data_list()` for anonymous S3 listings using
      `region = "eu-west-1"`.
- [ ] Implement `ont_open_data_fetch()` with explicit object-key selection,
      cache directory control, and no implicit bulk downloads.
- [ ] Add a named dataset helper for
      `s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/`.
- [ ] Add named constants or helper accessors for the selected example objects:
      `PAU85136_pass_279c9095_68316534_8289.pod5` and
      `PAU85136_fail_279c9095_68316534_0.pod5`.
- [ ] Record object key, size, last-modified timestamp, and local cache path in
      returned metadata.
- [ ] Add tests with mocked S3 listings and no network requirement.
- [ ] Document that downloaded POD5 files must not be committed.

## Slice 8: POD5 Verification And File Info

- [ ] Expose `pod5_verify(path)` through Rust-in-R.
- [ ] Expose `pod5_file_info(path)` through Rust-in-R.
- [ ] Return structured data frames/lists rather than TSV/JSON text.
- [ ] Map Rust path, format, schema, and integrity failures to R conditions.
- [ ] Add tests for valid POD5-like signatures, non-POD5 files, truncated files,
      and missing paths.

## Slice 9: POD5 Folder Info And Manifests

- [ ] Expose `pod5_folder_info(path)` for run-tree QC.
- [ ] Expose `pod5_manifest(path)` with a versioned schema.
- [ ] Expose `pod5_compare(left, right)` for operational handoff checks.
- [ ] Add tests for mixed flow cells, duplicate file names, failed files, and
      changed manifests.
- [ ] Document how these outputs feed synoptikon.

## Slice 10: Subdivision And Playback Planning

- [ ] Expose read-only `pod5_subdivide_plan()` in R.
- [ ] Expose playback planning only if synoptikon/mnemosyne needs simulated run
      arrivals for QC/review workflows.
- [ ] Keep write-capable POD5 subdivision out of R until read-only contracts are
      stable.
- [ ] Add deterministic tests for planning outputs.
- [ ] Define the demonstration subset policy for `pod5-tools`-derived POD5
      splits, including source provenance, subsetting parameters, output
      checksums, and artifact storage location.
- [ ] Add a small derived-POD5 workflow for examples once the read-only binding
      and provenance contracts are stable.

## Slice 11: BAM/Bamana Integration

- [ ] Audit `../bamana` public library APIs and identify only the QC/review
      command logic that must be promoted for in-process R use.
- [ ] Expose `bam_summary(path)` through Rust-in-R with stable R-native return
      shapes.
- [ ] Expose `bam_verify(path)`, `bam_validate(path)`, and `bam_check_eof(path)`
      through Rust-in-R.
- [ ] Expose BAM index, mapping, sorting, and tag checks needed for QC report
      cards.
- [ ] Map Bamana response envelopes and errors to R conditions without losing
      command semantics.
- [ ] Add small BAM fixtures and offline tests for summary/verify/validate
      wrappers.
- [ ] Define BAM QC report-card checks for mapping fraction, duplicate
      fraction, QC-fail fraction, MAPQ zero burden, stale/missing index,
      sorting mismatch, missing expected tags, EOF absence, validation
      findings, and provenance anomalies.
- [ ] Document which Bamana operations are read-only QC surfaces and which are
      transformations requiring explicit output paths.
- [ ] Explicitly document which Bamana capabilities are out of scope for
      floundeR because they do not serve nanopore QC/review directly.

## Slice 12: Synoptikon Export

- [ ] Define a versioned synoptikon QC payload.
- [ ] Implement `as_synoptikon_qc()` or `write_synoptikon_qc()`.
- [ ] Include sequencing-summary, flowcell, barcode, POD5 integrity, BAM QC,
      and library-preparation evidence sections.
- [ ] Add validation tests against the expected `../mnemosyne` schema.
- [ ] Add a vignette showing end-to-end handoff.

## Slice 12A: Porkchop Library-Preparation QC

- [ ] Audit `../porkchop` public library APIs and identify only the QC/review
      logic that should be promoted for in-process R use.
- [ ] Expose curated R wrappers for kit candidates, adapter/primer evidence,
      barcode evidence, and cDNA primer-pair evidence.
- [ ] Preserve Porkchop heuristic score terminology; do not expose scores as
      calibrated probabilities.
- [ ] Carry kit provenance, lifecycle status, support level, and validation
      limitations into floundeR return values.
- [ ] Add small FASTQ fixtures or mocked sequence records for offline tests.
- [ ] Define report-card checks for unexpected library chemistry, adapter
      burden, barcode ambiguity, cDNA partial/unclassified burden, and
      unsupported or partial kit support levels.
- [ ] Document which Porkchop capabilities are out of scope for floundeR because
      they are preprocessing/trimming features rather than QC/review evidence.

## Slice 13: Grammateus Reporting Integration

- [ ] Confirm the local `../grammateus` source checkout is current and clean
      before adding it as a Rust library dependency.
- [ ] Verify the Grammateus charter, SRS, ARDs, and rendering/trusted-report
      contracts in `../mnemosyne-docs` before changing `../grammateus`.
- [ ] Identify the Grammateus public library APIs needed for semantic report
      elements, rendering, manifests, provenance, branding, and trusted-report
      lifecycle metadata.
- [ ] Implement R-to-Grammateus plot artifact handoff for existing `ggplot2`
      plots, including deterministic SVG/PNG output, checksums, captions, alt
      text, methods notes, and provenance.
- [ ] Implement Grammateus-to-R controlled plot generation wrappers for
      semantic `ReportPlot` specs, using local `Rscript` or governed container
      execution as configured.
- [ ] Add in-process Rust bindings from `floundeR` to Grammateus report
      rendering; do not shell out to Grammateus CLI binaries.
- [ ] Define the prebuilt Grammateus runtime interface used by public floundeR
      builds when private Grammateus source is unavailable.
- [ ] Implement Grammateus runtime discovery helpers such as
      `grammateus_runtime_available()`, `grammateus_runtime_version()`,
      `grammateus_runtime_validate()`, and `grammateus_runtime_install()`.
- [ ] Ensure core QC APIs and package checks do not require Grammateus runtime
      assets.
- [ ] Define Grammateus semantic elements for run metadata, QC summary,
      flowcell density, yield over time, quality distribution, barcode balance,
      POD5 integrity, BAM alignment summaries, BAM validation/index/sort
      evidence, Porkchop kit/library-preparation evidence, report-card
      findings, methods, limitations, and appendices.
- [ ] Encode Mnemosyne Biosciences branding and styling through Grammateus
      templates/themes.
- [ ] Implement `qc_report()` or equivalent R API that produces HTML/PDF plus a
      report manifest/provenance payload.
- [ ] Add the first QC plot helper set for yield over time, read quality, read
      length, flowcell density, barcode balance, POD5 integrity, BAM mapping,
      BAM MAPQ, and BAM flag summaries.
- [ ] Add tests for report element schemas, required captions/alt text,
      provenance hashes, and stable report identifiers.
- [ ] Keep RMarkdown support only as a documented legacy/transitional path until
      Grammateus reports cover the required technical-report workflow.

## Slice 13A: GitHub And Future BioConductor Distribution

- [ ] Define GitHub release asset names and manifests for prebuilt Grammateus
      runtime bundles.
- [ ] Add checksum and signature verification for private Grammateus runtime
      artifacts.
- [ ] Add compatibility checks between floundeR versions and Grammateus runtime
      versions.
- [ ] Add CI separation: public floundeR checks without Grammateus, credentialed
      private checks with Grammateus runtime/report rendering.
- [ ] Document GitHub installation paths for open-source users and authorized
      Grammateus runtime users.
- [ ] Before any BioConductor submission, verify current policy for optional
      external binaries, private runtime assets, downloads, and system
      requirements.

## Slice 14: Documentation Refresh

- [ ] Rewrite README around the contemporary QC workflow.
- [ ] Add a "run health in one minute" example.
- [ ] Replace FAST5 vignette with POD5, BAM, and synoptikon QC vignettes.
- [ ] Use the ONT Zymo fecal POD5 dataset as the canonical opt-in real-data
      POD5 example.
- [ ] Use `pod5-tools`-derived subset/split POD5 artifacts for lightweight
      demonstrations where full source POD5 files are too large.
- [ ] Add a Grammateus-backed technical report vignette.
- [ ] Add a vignette showing both reporting directions: R plot dropped into a
      Grammateus report, and Grammateus semantic plot spec rendered by R.
- [ ] Update pkgdown navigation.
- [ ] Add plot examples for report-card findings.

## Slice 15: Release Preparation

- [ ] Bump version according to semantic versioning.
- [ ] Update `NEWS.md`.
- [ ] Run full tests.
- [ ] Run `R CMD check`.
- [ ] Build documentation.
- [ ] Commit and push changes after explicit prompt or as required by the active
      agent instructions.
- [ ] Defer tagging/release until explicit approval.
