# Bamana Integration Notes

This document records floundeR's curated Bamana integration decisions. Bamana
support must remain focused on nanopore QC, review, provenance, reporting, and
synoptikon handoff. It must not become a wholesale R mirror of Bamana's public
CLI or transformation surface.

## Audit Preflight

Checked on 2026-06-14:

- `../bamana/CHARTER.md`
- `../bamana/docs/project-charter.md`
- `../bamana/docs/json-output.md`
- `../bamana/docs/cli.md`
- `../bamana/spec/jsonschema/`
- `../bamana/spec/examples/`
- `../bamana/src/lib.rs`
- `../bamana/src/main.rs`
- `../bamana/src/commands/`

No Bamana-specific canonical files were found in `../mnemosyne-docs` during the
floundeR governance review, so the local Bamana charter, JSON contracts, and
public schema/examples remain the governing boundary until canonical
mnemosyne-docs material exists.

Local Bamana source state at audit time:

- crate: `bamana`
- crate version: `0.1.0`
- git commit: `d64c1d4`

## Public Contract Shape

Bamana is a JSON-first Rust toolkit. The public command envelope is governed by
`../bamana/spec/jsonschema/common/command_response.schema.json` and contains:

- `ok`
- `command`
- `path`
- `analysis_wall_seconds`
- `data`
- `error`

floundeR should preserve those response semantics internally, but R wrappers
should return idiomatic R data frames/lists and typed conditions rather than
JSON strings. The Rust-in-R layer should normalize Bamana's command responses
before converting to R objects so that failure payloads carrying partial data
are not lost.

The existing Bamana Rust library already exposes command modules such as
`bamana::commands::summary`, `verify`, `validate`, `check_eof`,
`check_index`, `check_map`, `check_sort`, `check_tag`, `checksum`, `header`,
`inspect_duplication`, and `forensic_inspect`. Their request and payload
structs are public in the command modules. Some commands return
`CommandResponse<T>` directly, while others return `Result<T, AppError>` and
are wrapped into `CommandResponse<T>` in `src/main.rs`. floundeR should add or
use a small Bamana library adapter that returns a consistent envelope for the
curated commands.

## First Curated Read-Only Surface

The first floundeR BAM surface should promote only read-only QC and review
contracts:

| floundeR target | Bamana source | Purpose |
| --- | --- | --- |
| `bam_summary(path)` | `commands::summary::run()` | Operational alignment overview: counts, mapping fractions, MAPQ, optional flags, index-derived evidence, and bounded/full-scan evidence scope. |
| `bam_verify(path)` | `commands::verify::run()` | Header-level BAM/BGZF verification only; this does not prove EOF completeness or alignment-record validity. |
| `bam_validate(path)` | `commands::validate::run()` | Structural validation findings with explicit header-only, bounded-record, or full validation mode. |
| `bam_check_eof(path)` | `commands::check_eof::run()` | Canonical BGZF EOF-marker evidence only. |
| `bam_check_index(path)` | `commands::check_index::run()` | Adjacent index discovery, BAI/CSI/GZI support level, usability, stale state, and candidate sidecars. |
| `bam_check_map(path)` | `commands::check_map::run()` | Mapping-state evidence from usable index metadata or scanner fallback, including region-scoped evidence where requested later. |
| `bam_check_sort(path)` | `commands::check_sort::run()` | Declared versus observed sort evidence with bounded/strict evidence scope. |
| `bam_check_tag(path, tag)` | `commands::check_tag::run()` | Auxiliary tag presence evidence with bounded/full-scan semantics. |
| `bam_checksum(path)` | `commands::checksum::run()` | Explicit checksum domains for provenance and handoff comparisons. |
| `bam_header(path)` | `commands::header::run()` | Header/reference/read-group/program metadata and non-fatal reference diagnostics. |
| `bam_forensic_inspect(path)` | `commands::forensic_inspect::run()` | Provenance and operational anomaly evidence, not fraud accusation. |
| `bam_inspect_duplication(path)` | `commands::inspect_duplication::run()` | Collection-duplication evidence, distinct from PCR duplicate marking and BAM duplicate flags. |

`bam_summary()`, `bam_verify()`, `bam_validate()`, and `bam_check_eof()` are the
minimum first binding set because they satisfy the next Slice 11 TODO items and
provide the core report-card inputs for aligned-read QC.

`bam_summary()` is implemented in floundeR 0.6.0 as the first binding in this
set. It links Bamana as an in-process Rust library dependency pinned to
`d64c1d4cba524c4ef24bc45c6b9b721881b5009c`, calls
`commands::summary::run()`, and returns a named R list of schema-stable data
frames for status, evidence, header, counts, fractions, MAPQ, mapping,
anomalies, flag categories, references, index-derived evidence, and optional
MAPQ histograms. It deliberately leaves region-scoped summaries and the wider
verify/validate/index/tag surfaces for subsequent Slice 11 bindings.

## R Return-Shape Guidance

R wrappers should expose stable, schema-versioned objects rather than
unstructured nested JSON. The first binding set should use:

- a one-row command/status data frame with schema version, command, path, `ok`,
  analysis seconds, evidence mode, and semantic note;
- focused tidy child data frames for counts, fractions, MAPQ, validation
  findings, index candidates, tag results, and checksum domains;
- a typed R condition family for Bamana errors, carrying Bamana `error.code`,
  `message`, `detail`, and `hint`;
- explicit columns that distinguish bounded evidence from full-file evidence,
  because bounded scans must not be interpreted as complete-file claims.

The R layer should not present Bamana validation, mapping, checksum, or forensic
results as biological truth. It should carry Bamana's own semantic notes into
report elements and synoptikon payloads.

## Transformation Boundary

The following Bamana capabilities are useful but should not be exposed in
floundeR's first BAM QC API because they write outputs or alter records:

- `sort`
- `merge`
- `explode`
- `select_region`
- `consume`
- `filter`
- `subsample`
- `fastq`
- `unmap`
- `reheader`
- `annotate_rg`
- `deduplicate`
- `index`

These operations may be added later only when they directly support QC/review,
require explicit output paths, and return provenance-rich reports. They should
remain distinct from read-only QC wrappers and must preserve Bamana's command
semantics, output-safety rules, and JSON contracts.

## Report And Synoptikon Use

Bamana-derived floundeR report sections should include:

- BAM identity, shallow verification, EOF completeness, and validation state;
- mapping summary, MAPQ zero burden, duplicate and QC-fail fractions, and
  primary/secondary/supplementary composition;
- index state, stale/missing index findings, and scan/index evidence source;
- sort-order agreement;
- expected tag evidence such as read group tags where requested;
- checksum/provenance domains for handoff comparison;
- forensic and collection-duplication findings when enabled.

These outputs should feed future Grammateus report elements and the planned
synoptikon QC payload without requiring private Grammateus runtime assets.
