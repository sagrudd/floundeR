# Rust-In-R Integration Decision

This document records the Slice 6 architecture decision for compiled
Rust support in `floundeR`.

## Decision

`floundeR` will use one embedded Rust extension crate owned by this
repository, with R bindings exposed through `extendr`.

Development-time Rust dependencies will be path dependencies on adjacent
repositories:

- `../pod5-tools` for POD5 discovery, verification, file information,
  folder-level summaries, manifests, comparison, subdivision planning,
  and demonstration-subset provenance.
- `../bamana` for BAM/BGZF/FASTQ summary, verification, validation,
  EOF/index evidence, mapping/sort/tag checks, and forensic/provenance
  signals.
- `../porkchop` for kit-registry metadata, adapter/primer/barcode
  evidence, cDNA primer-pair evidence, and library-preparation mismatch
  signals.
- `../grammateus` only for authorized private development builds of
  report element, rendering, manifest, provenance, branding, and
  trusted-report lifecycle behavior.

The public package must remain installable and checkable without private
Grammateus source. Public distribution will use optional prebuilt
Grammateus runtime/reporting assets as described in `DISTRIBUTION.md`,
not a required source dependency.

## Rejected Shape

Separate R support packages are not the initial target. They would add
package release coordination before the core floundeR QC contract is
stable. A split package can be reconsidered only if compiled dependency
size, Bioconductor policy, or platform binary distribution makes the
embedded crate untenable.

CLI wrappers are also rejected for normal package APIs. Shelling out to
`pod5-tools`, `bamana`, `porkchop`, or Grammateus may be useful for
development diagnostics, but user-facing floundeR APIs must call Rust
library functions in-process and return R-native data frames, lists, and
typed R conditions.

## Curation Rule

Rust-backed APIs may be added to floundeR only when they directly
support at least one of these QC and review needs:

- nanopore run QC;
- raw-data or alignment integrity review;
- provenance, manifests, comparison, or operational handoff;
- report-card or technical-report evidence;
- Synoptikon/Mneion trusted lifecycle integration;
- small demonstration datasets or derived-fixture provenance.

Rust-backed APIs should not be added merely because an upstream engine
supports the operation. In particular:

- POD5 APIs must stay focused on discovery, integrity, metadata,
  manifests, comparison, subdivision planning, and demonstration
  provenance.
- BAM/Bamana APIs must stay focused on read-only QC evidence,
  validation, mapping/index/sort/tag/EOF state, summaries, and explicit
  provenance.
- Porkchop APIs must stay focused on library-preparation evidence and
  must preserve heuristic score terminology rather than presenting
  scores as calibrated probabilities.
- Grammateus APIs must stay focused on semantic report elements,
  rendering, provenance, manifests, and trusted-report lifecycle
  payloads. floundeR must not build a parallel report trust stack.

Transformations that write output files, such as filtering, trimming,
subsampling, sorting, merging, reheadering, POD5 subdivision, or BAM
annotation, require explicit output paths and provenance reporting. They
should not be introduced before the corresponding read-only QC contracts
are stable.

## R Boundary

The R layer owns analyst ergonomics:

- exported functions;
- input validation;
- conversion to tibbles/lists;
- stable schema-version labels;
- condition classes and actionable error messages;
- test skipping when compiled Rust support is unavailable.

The Rust layer owns performance-sensitive parsing, scanning, validation,
manifest generation, and renderer calls. Rust structs must be converted
into R-native shapes before being returned to users.

## Dependency Boundary

The embedded crate should depend on open engines by path during local
development and by released git/crate revisions when public distribution
is prepared. The release manifest should record tested versions of
`pod5-tools`, `bamana`, and `porkchop`.

Grammateus remains private. The public package should expose runtime
discovery and validation helpers for prebuilt Grammateus assets, but
core QC APIs and package checks must not require those assets.

## Contract Checks Before Cross-Repository Changes

Before changing an adjacent repository, inspect and preserve:

- `../pod5-tools/AGENTS.md`, `../pod5-tools/README.md`,
  `../pod5-tools/roadmap.md`, and `../pod5-tools/todo.md`;
- `../bamana/CHARTER.md`, `../bamana/docs/project-charter.md`, and
  Bamana public JSON/schema documentation;
- `../porkchop/AGENTS.md`, `../porkchop/ROADMAP.md`,
  `../porkchop/README.md`, and the relevant Sphinx provenance/output/
  validation documents;
- the Grammateus and Mneion contracts listed in `GOVERNANCE.md`.

If an intended floundeR integration conflicts with those documents,
update the plan rather than silently widening the integration surface.

## Current Cross-Repository Preflight

On 2026-06-13, the Slice 6 scaffold preflight inspected the governing
documents listed above for pod5-tools, Bamana, Porkchop, Grammateus, and
Mneion/Synoptikon. The first `extendr` scaffold target is a local
capability/version function and does not require changes to
`../pod5-tools`, `../bamana`, or `../porkchop`.

Therefore, the next scaffold slice may proceed without adjacent
repository edits. If a later functional binding requires an adjacent
change, that binding slice must re-read the relevant governing
documents, make the smallest library API change in that repository, run
that repository’s relevant checks, and preserve its public contract and
terminology.

## First Scaffold Target

The first compiled scaffold should be deliberately minimal:

1.  Add the `extendr` build skeleton.
2.  Expose a small capability/version function from Rust.
3.  Add R wrappers that skip cleanly unless compiled support is
    available.
4.  Keep all existing pure-R QC contracts operational without compiled
    Rust.

Only after that scaffold is checked should
[`pod5_find()`](https://sagrudd.github.io/floundeR/reference/pod5_find.md)
become the first functional Rust-backed QC API.

## Scaffold Layout

The initial scaffold lives under `src/` so it participates in normal R
source package builds:

- `src/rust/Cargo.toml` defines the private `flounder-extendr` static
  library crate using `extendr-api` plus curated open Rust QC
  dependencies.
- `src/rust/src/lib.rs` defines registered R-callable Rust entry points.
- `src/Makevars` and `src/Makevars.win` build the Rust static library
  before the R shared object is linked.
- `src/init.c` provides a conservative native registration entry point.

The extension links curated POD5 and BAM/BGZF/FASTQ dependencies through
the open `pod5-tools` and `bamana` Git repositories. Porkchop-backed R
wrappers are present, but the public default build keeps Porkchop
unlinked until the Porkchop repository is public; local sibling-checkout
builds use `src/rust-porkchop/Cargo.toml` and the `porkchop-integration`
feature. Private Grammateus should be introduced only through the
optional runtime model documented for reporting.

## Minimal Callable Function

The first Rust-to-R proof point was an internal capability function
exposed as `.flounder_rust_capabilities()`. The current Rust entry
points are registered through `src/init.c`; capability metadata reports
linked status for `pod5-tools` and `bamana`, feature-dependent status
for Porkchop, and explicit `not_linked` status for Grammateus. The
internal helpers remain in place so public functions can enforce a
consistent R error and availability policy.

## R Wrapper And Error Policy

The public R wrapper surface now starts with
[`flounder_rust_capabilities()`](https://sagrudd.github.io/floundeR/reference/flounder_rust_capabilities.md),
[`flounder_rust_available()`](https://sagrudd.github.io/floundeR/reference/flounder_rust_capabilities.md),
and
[`skip_if_no_flounder_rust()`](https://sagrudd.github.io/floundeR/reference/flounder_rust_capabilities.md).
Package code should call these wrappers rather than calling
[`.Call()`](https://rdrr.io/r/base/CallExternal.html) directly. Missing
compiled support is represented as a typed `floundeR_rust_unavailable`
condition when a Rust-backed feature is required, or as capability
metadata with `compiled_support = FALSE` when callers are only probing
availability.

Future POD5, Bamana, Porkchop, and Grammateus bindings should follow the
same pattern:

- public R functions validate inputs and return R-native objects;
- internal native wrappers convert Rust failures into typed R
  conditions;
- tests that require optional compiled or external runtime support call
  [`skip_if_no_flounder_rust()`](https://sagrudd.github.io/floundeR/reference/flounder_rust_capabilities.md)
  or a more specific helper built on the same capability metadata.

## Source Install Requirements

Source installs require Cargo and rustc because `src/Makevars` builds
the embedded Rust static library before linking the R package shared
object. The current minimum Rust version is `1.85`. Detailed macOS,
Linux, Windows, Docker, and CI setup notes live in
`DEVELOPMENT.md#source-install-requirements`.

The public package must not require private Grammateus source or runtime
assets to install or check. Public CI should verify the open-source Rust
scaffold and core QC package only; private Grammateus report rendering
belongs in a separate credentialed check path.

## POD5 Discovery Preflight

Slice 7 inspected `../pod5-tools` on 2026-06-13 before adding the first
functional POD5 binding. The existing `pod5-tools` library already
exposes the curated discovery API floundeR needs:

- `pod5_tools::find_pod5_directories(root: PathBuf)`;
- `pod5_tools::Pod5DirectoryRecord` with path, POD5 file count, byte
  total, oldest modification time, and newest modification time.

No adjacent `pod5-tools` code promotion is required for the initial
[`pod5_find()`](https://sagrudd.github.io/floundeR/reference/pod5_find.md)
binding. The floundeR wrapper should call that library function
in-process and convert records to a tibble/data frame, rather than
invoking the `pod5-tools find` CLI or parsing TSV/JSON text.

The toolchain prerequisite is now aligned: floundeR requires Rust `1.85`
or newer and the embedded crate uses Rust edition 2024, matching the
floor needed to depend on the current `../pod5-tools` crate.

For public source-package installs, the embedded crate pins `pod5-tools`
to the open GitHub repository at the commit matching the inspected
adjacent checkout. That avoids brittle release-tarball failures caused
by local sibling path dependencies while preserving the same in-process
Rust library boundary.

## Compiled-Code Check Warning

The first `pod5-tools` binding necessarily moved the embedded crate from
the minimal `no_std` proof point to a standard Rust library stack
because POD5 discovery walks the filesystem, owns path buffers, and
reports modification timestamps through `pod5-tools`. In the Docker
release-tarball check, R reports a compiled-code warning for `_exit`,
`abort`, and `exit` symbols in `floundeR.so` after linking the Rust
static library.

floundeR does not intentionally call process-terminating functions from
its public R APIs. Rust failures must continue to be mapped into typed R
conditions or data-frame/list responses at the R boundary. The warning
is therefore tracked as release-readiness debt for the Rust integration
stack rather than as a blocker for the current GitHub-source reboot
line.

Before Bioconductor/CRAN-facing release work, this must be revisited.
Acceptable outcomes are:

- a supported Rust/R build strategy that removes the warning while
  preserving in-process library calls;
- explicit policy acceptance recorded in the release checklist;
- a split binary/runtime distribution design if package-policy
  constraints make embedded Rust static linking untenable.

Do not replace Rust library calls with CLI subprocess wrappers merely to
avoid this warning; that would violate the core architecture decision.
