# Rust-In-R Integration Decision

This document records the Slice 6 architecture decision for compiled Rust
support in `floundeR`.

## Decision

`floundeR` will use one embedded Rust extension crate owned by this repository,
with R bindings exposed through `extendr`.

Development-time Rust dependencies will be path dependencies on adjacent
repositories:

- `../pod5-tools` for POD5 discovery, verification, file information,
  folder-level summaries, manifests, comparison, subdivision planning, and
  demonstration-subset provenance.
- `../bamana` for BAM/BGZF/FASTQ summary, verification, validation, EOF/index
  evidence, mapping/sort/tag checks, and forensic/provenance signals.
- `../porkchop` for kit-registry metadata, adapter/primer/barcode evidence,
  cDNA primer-pair evidence, and library-preparation mismatch signals.
- `../grammateus` only for authorized private development builds of report
  element, rendering, manifest, provenance, branding, and trusted-report
  lifecycle behavior.

The public package must remain installable and checkable without private
Grammateus source. Public distribution will use optional prebuilt Grammateus
runtime/reporting assets as described in `DISTRIBUTION.md`, not a required
source dependency.

## Rejected Shape

Separate R support packages are not the initial target. They would add package
release coordination before the core floundeR QC contract is stable. A split
package can be reconsidered only if compiled dependency size, Bioconductor
policy, or platform binary distribution makes the embedded crate untenable.

CLI wrappers are also rejected for normal package APIs. Shelling out to
`pod5-tools`, `bamana`, `porkchop`, or Grammateus may be useful for development
diagnostics, but user-facing floundeR APIs must call Rust library functions
in-process and return R-native data frames, lists, and typed R conditions.

## Curation Rule

Rust-backed APIs may be added to floundeR only when they directly support at
least one of these QC and review needs:

- nanopore run QC;
- raw-data or alignment integrity review;
- provenance, manifests, comparison, or operational handoff;
- report-card or technical-report evidence;
- Synoptikon/Mneion trusted lifecycle integration;
- small demonstration datasets or derived-fixture provenance.

Rust-backed APIs should not be added merely because an upstream engine supports
the operation. In particular:

- POD5 APIs must stay focused on discovery, integrity, metadata, manifests,
  comparison, subdivision planning, and demonstration provenance.
- BAM/Bamana APIs must stay focused on read-only QC evidence, validation,
  mapping/index/sort/tag/EOF state, summaries, and explicit provenance.
- Porkchop APIs must stay focused on library-preparation evidence and must
  preserve heuristic score terminology rather than presenting scores as
  calibrated probabilities.
- Grammateus APIs must stay focused on semantic report elements, rendering,
  provenance, manifests, and trusted-report lifecycle payloads. floundeR must
  not build a parallel report trust stack.

Transformations that write output files, such as filtering, trimming,
subsampling, sorting, merging, reheadering, POD5 subdivision, or BAM annotation,
require explicit output paths and provenance reporting. They should not be
introduced before the corresponding read-only QC contracts are stable.

## R Boundary

The R layer owns analyst ergonomics:

- exported functions;
- input validation;
- conversion to tibbles/lists;
- stable schema-version labels;
- condition classes and actionable error messages;
- test skipping when compiled Rust support is unavailable.

The Rust layer owns performance-sensitive parsing, scanning, validation,
manifest generation, and renderer calls. Rust structs must be converted into
R-native shapes before being returned to users.

## Dependency Boundary

The embedded crate should depend on open engines by path during local
development and by released git/crate revisions when public distribution is
prepared. The release manifest should record tested versions of `pod5-tools`,
`bamana`, and `porkchop`.

Grammateus remains private. The public package should expose runtime discovery
and validation helpers for prebuilt Grammateus assets, but core QC APIs and
package checks must not require those assets.

## Contract Checks Before Cross-Repository Changes

Before changing an adjacent repository, inspect and preserve:

- `../pod5-tools/AGENTS.md`, `../pod5-tools/README.md`,
  `../pod5-tools/roadmap.md`, and `../pod5-tools/todo.md`;
- `../bamana/CHARTER.md`, `../bamana/docs/project-charter.md`, and Bamana
  public JSON/schema documentation;
- `../porkchop/AGENTS.md`, `../porkchop/ROADMAP.md`,
  `../porkchop/README.md`, and the relevant Sphinx provenance/output/
  validation documents;
- the Grammateus and Mneion contracts listed in `GOVERNANCE.md`.

If an intended floundeR integration conflicts with those documents, update the
plan rather than silently widening the integration surface.

## First Scaffold Target

The first compiled scaffold should be deliberately minimal:

1. Add the `extendr` build skeleton.
2. Expose a small capability/version function from Rust.
3. Add R wrappers that skip cleanly unless compiled support is available.
4. Keep all existing pure-R QC contracts operational without compiled Rust.

Only after that scaffold is checked should `pod5_find()` become the first
functional Rust-backed QC API.
