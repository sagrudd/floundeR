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

`../pod5-tools` is currently a Rust edition 2024 crate. floundeR still declares
Rust `1.71` as its minimum. Before adding `../pod5-tools` as a path dependency,
align floundeR's Rust toolchain floor with the upstream crate and update the
source-install documentation, Docker/CI expectations, and package metadata.

## Scope Guard

For the first POD5 binding, floundeR should expose only `pod5_find(path)`.
Verification, file information, folder summaries, manifests, comparison, and
subdivision planning remain later slices. Write-capable POD5 subdivision stays
out of the R API until read-only contracts are stable.
