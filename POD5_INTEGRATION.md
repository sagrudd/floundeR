# POD5 Integration Notes

This document records floundeR's curated POD5 integration decisions. POD5
support must remain focused on nanopore QC, review, provenance, reporting,
demonstration datasets, and synoptikon handoff. It must not become a wholesale
R mirror of `../pod5-tools` or the official ONT POD5 toolkit.

## Discovery Preflight

Checked on 2026-06-13:

- `../pod5-tools/AGENTS.md`
- `../pod5-tools/README.md`
- `../pod5-tools/roadmap.md`
- `../pod5-tools/todo.md`
- `../pod5-tools/src/lib.rs`

The discovery logic needed for floundeR is already available as a pure Rust
library API in `../pod5-tools`:

- `pod5_tools::find_pod5_directories(root: PathBuf)`
- `pod5_tools::Pod5DirectoryRecord`

The record contains the exact Slice 7 fields required for `pod5_find()`:
directory path, POD5 file count, total bytes, oldest modified timestamp, and
newest modified timestamp. The focused upstream `find` tests pass with:

```sh
cd ../pod5-tools
cargo test find --lib
```

No `pod5-tools` code promotion is needed for the first floundeR discovery
binding. The correct floundeR implementation is an in-process Rust call to the
library function, followed by conversion to an R data frame or tibble.

## Dependency Prerequisite

`../pod5-tools` is currently a Rust edition 2024 crate. floundeR now requires
Rust `1.85` or newer and the embedded crate also uses Rust edition 2024, so the
toolchain floor is aligned for the first `../pod5-tools` path dependency.

## Verification And File Info

Slice 8 added in-process Rust bindings for:

- `pod5_tools::verify_pod5_file(path)`
- `pod5_tools::read_pod5_file_info(reader, path)`
- `pod5_tools::FilesystemPod5MetadataReader`

The R functions `pod5_verify()` and `pod5_file_info()` return data frames and
typed R conditions. The current `pod5-tools` backend performs fast extension and
signature verification and reports deeper parser-backed checks as unavailable or
`not_checked` until the concrete POD5 parser backend is connected.

## Folder Info, Manifests, And Compare

Slice 9 added in-process Rust bindings for:

- `pod5_tools::folder_info(path, reader)`
- `pod5_tools::manifest_from_path(path)`
- `pod5_tools::compare_inputs(left, right)`

The R functions `pod5_folder_info()`, `pod5_manifest()`, and `pod5_compare()`
are read-only QC and handoff surfaces. They intentionally expose a curated
subset of `pod5-tools`: run-tree folder summaries, versioned collection
manifests, and operational drift checks. They do not expose POD5 rewriting,
subdivision writing, playback emission, or general ONT tool replacement
behavior.

Mixed flow-cell and mixed sequencing-kit aggregation is covered in
`pod5-tools` by a path-aware mocked `Pod5MetadataReader` test. The public
floundeR R binding currently uses `FilesystemPod5MetadataReader`, which cannot
recover real flow-cell or kit metadata from POD5 internals. Public R-level
mixed-flow-cell tests should be added when the parser-backed reader is exposed
through the same binding; until then, `pod5_folder_info()` reports unavailable
metadata explicitly in `flow_cell_ids`, `sequencing_kits`, and `warnings`.

## Synoptikon Handoff

The POD5 collection outputs are designed to feed Synoptikon/Mnemosyne QC
payloads as follows:

- `pod5_folder_info()` supplies the run-tree raw-data integrity summary:
  POD5 file count, byte total, available read count, duplicate file names,
  verification failure count, integrity status, and warnings.
- `pod5_manifest()` supplies the versioned POD5 file inventory for provenance:
  manifest schema version, source path, relative paths, file sizes,
  verification status, and failed-check counts.
- `pod5_compare()` supplies operational handoff drift evidence:
  files missing on either side, changed file sizes, changed verification
  status, and an explicit `match` row when collections agree.

These outputs should be included in the future `as_synoptikon_qc()` payload
under a POD5/raw-data integrity section. They should also be eligible report
inputs for Grammateus POD5 integrity plots and provenance tables.

## Subdivision Planning

Slice 10 added an in-process Rust binding for:

- `pod5_tools::subdivide_plan_from_path(path, strategy, files_per_chunk,
  seconds_per_chunk, reads_per_chunk)`

The R function `pod5_subdivide_plan()` is read-only. It returns a deterministic
plan table describing how a local POD5 file, folder, run tree, or manifest
could be grouped by file count, sample label, elapsed time, or read count. It
does not write POD5 files, copy files, create output directories, or expose the
`pod5-tools` write-capable subdivision API.

The current filesystem-backed POD5 metadata reader cannot recover acquisition
timestamps or read counts, so elapsed-time and read-count strategies return a
single placeholder chunk with an explicit warning. Those strategies should
become fully data-aware when the parser-backed POD5 reader is exposed through
the binding.

Playback planning remains out of floundeR's public API for now. It is useful in
`pod5-tools` for workflow development, but floundeR should only expose it after
synoptikon or `../mnemosyne` has a concrete QC/review requirement for simulated
run arrivals.

## Scope Guard

floundeR should keep exposing only read-only POD5 functions that directly
support QC, review, provenance, reporting, demonstration datasets, or synoptikon
handoff. Write-capable POD5 subdivision stays out of the R API until read-only
planning and provenance contracts are stable.
