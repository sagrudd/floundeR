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

`bam_summary()`, `bam_verify()`, `bam_validate()`, `bam_check_eof()`,
`bam_check_index()`, `bam_check_map()`, `bam_check_sort()`, and
`bam_check_tag()` are the initial read-only binding set because they provide
the core report-card inputs for aligned-read QC without exposing Bamana's
record-writing or transformation operations.

`bam_summary()` is implemented in floundeR 0.6.0 as the first binding in this
set. It links Bamana as an in-process Rust library dependency pinned to
`d64c1d4cba524c4ef24bc45c6b9b721881b5009c`, calls
`commands::summary::run()`, and returns a named R list of schema-stable data
frames for status, evidence, header, counts, fractions, MAPQ, mapping,
anomalies, flag categories, references, index-derived evidence, and optional
MAPQ histograms. It deliberately leaves region-scoped summaries and the wider
verify/validate/index/tag surfaces for subsequent Slice 11 bindings.

`bam_verify()`, `bam_validate()`, and `bam_check_eof()` are implemented in
floundeR 0.7.0. `bam_verify()` exposes Bamana's header-level verification as a
one-row evidence table and keeps EOF/body validation explicitly out of scope.
`bam_validate()` returns status, summary, findings, and error-metadata tables;
validation failures with findings are returned as QC evidence rather than
ordinary R exceptions. `bam_check_eof()` reports canonical BGZF EOF-marker
evidence as a one-row table, including missing EOF as `complete = FALSE`
evidence with Bamana error metadata.

`bam_check_index()`, `bam_check_map()`, `bam_check_sort()`, and
`bam_check_tag()` are implemented in floundeR 0.8.0. They expose Bamana's
read-only report-card surfaces for index availability/freshness, mapped-read
evidence, declared-versus-observed sort evidence, and aux-tag presence. The R
wrappers intentionally keep the argument surface narrow and return stable
tables for report assembly rather than mirroring Bamana's full command-line
surface.

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

Envelope semantics are part of the floundeR contract:

- Bamana responses with a structured payload are returned to R even when the
  Bamana command envelope has `ok = FALSE`. The returned R object preserves the
  command status in `status$ok` or `ok`, and carries any Bamana error metadata
  in an `error` table or `error_*` columns. This applies to QC evidence such as
  missing required indexes, validation findings, and missing BGZF EOF evidence.
- Bamana failures without a structured payload become typed R conditions in the
  `floundeR_bam_error` family. The condition preserves Bamana `code`,
  `message`, `detail`, and `hint` fields so callers can branch on stable
  metadata rather than parsing console text.
- R-side argument validation uses the same condition family with
  `floundeR_bam_argument_error`, so invalid caller input is separated from file,
  format, and index evidence.

`bam_qc_report_card()` is implemented in floundeR 0.9.0 as the first
Bamana-derived report-card contract. It does not call Bamana directly; instead
it consumes the stable R-native evidence tables returned by the curated BAM
wrappers. The first check set covers mapped-read fraction, duplicate-record
fraction, QC-fail fraction, MAPQ-zero burden, missing or unusable index, stale
index, sorting mismatch, missing expected aux tags, missing BGZF EOF marker,
validation finding count, and provenance anomaly count. Categorical checks are
encoded as 0/1 indicators so they fit the same pass/warn/fail schema used by
sequencing-summary report cards.

## Transformation Boundary

floundeR must separate read-only evidence from commands that write new files,
alter records, or normalize adjacent formats. The R package may call read-only
surfaces directly as QC evidence. Transformation surfaces may be added only as
explicit, output-path-bearing workflows with provenance-rich return objects.

| Bamana operation | floundeR posture | Boundary |
| --- | --- | --- |
| `identify` | read-only candidate | File identity evidence only; useful for provenance but not yet required by the first BAM card. |
| `verify` | read-only implemented | Header-level BAM/BGZF verification; does not prove EOF completeness or body validity. |
| `check_eof` | read-only implemented | BGZF EOF-marker evidence only. |
| `check_index` | read-only implemented | Adjacent sidecar discovery, support level, usability, and timestamp-staleness evidence. |
| `check_map` | read-only implemented | Mapping-state evidence from usable index metadata or bounded/full scan. Region-scoped evidence must remain distinct from whole-file evidence. |
| `check_sort` | read-only implemented | Declared versus observed sort evidence; not record rewriting. |
| `check_tag` | read-only implemented | Aux-tag presence/type evidence; bounded non-observation is not full-file absence unless full scan is reported. |
| `summary` | read-only implemented | Operational alignment overview; bounded summaries must not be interpreted as complete-file claims. |
| `validate` | read-only implemented | Structural and consistency findings over declared scope; not biological correctness. |
| `checksum` | read-only candidate | Handoff/provenance checksum domains; acceptable when surfaced as evidence and not as a write step. |
| `header` | read-only candidate | Header/reference/read-group/program/comment evidence; header diagnostics are not body validation. |
| `forensic_inspect` | read-only candidate | Provenance and operational anomaly evidence; must not be presented as fraud accusation. |
| `inspect_duplication` | read-only candidate | Collection-duplication evidence; distinct from PCR duplicate marking and BAM duplicate flags. |
| `benchmark` | out of normal floundeR API | Engineering/performance evidence, not routine biological or operational QC output. |
| `sort` | transformation | Writes a reordered BAM; requires explicit output path, input/output provenance, and index invalidation or regeneration notes. |
| `merge` | transformation | Writes a merged BAM; requires explicit output path, input ordering/compatibility provenance, and checksum reporting. |
| `explode` | transformation | Writes shards; requires explicit output directory, shard manifest, boundary policy, and checksums. |
| `select_region` | transformation | Writes selected records; requires explicit region semantics, duplicate-region policy, output path, header preservation, and provenance. |
| `consume` | transformation | Normalizes BAM/SAM/CRAM/FASTQ-like inputs into BAM; requires explicit output path and source-format/reference-policy provenance. |
| `filter` | transformation | Writes a retained-record BAM; requires explicit predicates, output path, and dropped/retained evidence. |
| `subsample` | transformation | Writes sampled BAM/FASTQ outputs; requires explicit seed/deterministic identity policy, output path, and sampling guarantees. |
| `fastq` | transformation | Exports FASTQ.GZ from BAM; loses alignment/header semantics by design and requires explicit output path. |
| `unmap` | transformation | Writes an unmapped BAM; strips mapping state and requires explicit output path and metadata-loss note. |
| `reheader` | transformation | Header-only mutation; must remain distinct from record-level `RG:Z` annotation and require explicit output path. |
| `annotate_rg` | transformation | Record-level read-group tag mutation; must remain distinct from header-only metadata mutation and require explicit output path. |
| `deduplicate` | transformation | Conservative collection-duplication remediation; distinct from PCR duplicate marking and requires dry-run/provenance-first reporting. |
| `index` | sidecar writer | Writes BAI/FASTQ.GZI sidecars; acceptable later only as an explicit write-capable helper with force/overwrite semantics and sidecar provenance. |

Write-capable operations must not be hidden behind report-card helpers or
implicit convenience calls. If floundeR adds any transformation helper later,
it must:

- require an explicit output file or output directory;
- return a structured provenance manifest with input/output paths, checksums,
  command parameters, Bamana version, floundeR version, and write status;
- preserve Bamana's distinction between dry-run planning and completed writes;
- report index invalidation, sidecar creation, or deferred indexing behavior;
- keep transformation outputs separate from read-only QC evidence in reports
  and synoptikon payloads.

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
