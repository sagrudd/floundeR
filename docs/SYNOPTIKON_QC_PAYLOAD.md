# Synoptikon QC Payload

This document defines the first floundeR-to-Synoptikon QC payload
contract. It is intentionally a contract document, not the exporter
implementation. The future
[`as_synoptikon_qc()`](https://sagrudd.github.io/floundeR/reference/as_synoptikon_qc.md)
or
[`write_synoptikon_qc()`](https://sagrudd.github.io/floundeR/reference/as_synoptikon_qc.md)
helper must emit payloads compatible with the installed JSON schema at:

``` text
inst/schema/synoptikon-qc-payload-v1.schema.json
```

## Contract Identity

- schema id:
  `https://mnemosyne.bio/schemas/flounder/synoptikon-qc-payload-v1.schema.json`
- schema version: `flounder.synoptikon_qc_payload.v1`
- producer: `floundeR`
- audience: Synoptikon/Mneion QC ingestion, review, governed reporting,
  and downstream audit evidence.

This payload is not a trusted report manifest and does not replace
Grammateus trusted-report lifecycle contracts. Grammateus remains
responsible for report rendering, report manifests, report provenance
tables, and signature lifecycle. The floundeR payload provides versioned
QC evidence that Synoptikon and Grammateus can consume.

## Top-Level Shape

The v1 payload is a JSON object with these required top-level fields:

- `schema_version`: exactly `flounder.synoptikon_qc_payload.v1`.
- `payload_id`: caller-supplied stable identifier for this QC handoff.
- `generated_at_utc`: ISO-8601 UTC production timestamp.
- `producer`: floundeR and runtime version metadata.
- `run_identity`: run, sample, project, flow-cell, kit, and instrument
  identifiers where known.
- `governance_context`: optional Synoptikon/Mneion identifiers such as
  tenant, governance domain, project, work request, work attempt, and
  report ids.
- `input_provenance`: array of local, object-store, or derived input
  records.
- `qc_sections`: sectioned QC evidence.
- `report_cards`: one or more pass/warn/fail report-card rows.
- `limitations`: explicit caveats, missing evidence, bounded scan notes,
  and parser/runtime limitations.
- `handoff`: intended consumer and compatibility metadata.

## QC Sections

`qc_sections` contains named section objects. Every section has:

- `present`: whether evidence for that section is present.
- `schema_version`: the native floundeR schema version for the evidence,
  or `null` when the section is not present.
- `status`: one of `pass`, `warn`, `fail`, `not_available`, or
  `not_checked`.
- `summary`: compact scalar values suitable for list pages or triage.
- `tables`: named arrays of row objects using existing floundeR table
  shapes.
- `artifacts`: optional references to external table, plot, manifest, or
  report artifacts with checksums.

The v1 named sections are:

- `sequencing_summary`: run-level read, base, quality, duration, and
  pass/fail evidence from
  [`qc_run_summary()`](https://sagrudd.github.io/floundeR/reference/qc_run_summary.md)
  and the tidy QC helpers.
- `flowcell`: channel-density and spatial-yield evidence.
- `barcode`: barcode composition and balance evidence.
- `pod5`: POD5 discovery, verification, folder-info, manifest,
  comparison, and demonstration-subset provenance evidence.
- `bam`: Bamana-derived summary, verification, validation, EOF, index,
  mapping, sorting, tag, and BAM report-card evidence.
- `library_preparation`: reserved for Porkchop-derived kit, adapter,
  primer, barcode, cDNA, score-terminology, and provenance evidence.
  Until Slice 12A lands, this section should normally be
  `present = false` and `status = "not_available"`.

## Provenance Rules

Input provenance records must not expose private storage credentials or
raw tenant storage topology. Use governed identifiers when a
Synoptikon/Mneion context exists, and include hashes, sizes, object
keys, or local paths only when they are policy-appropriate for the
execution context.

ONT open-data examples should record:

- bucket: `ont-open-data`
- region: `eu-west-1`
- object key
- object size
- last-modified timestamp
- local cache path only when the payload is local/operator-facing
- checksum or manifest hash when available

Derived POD5 demonstrations must record the source object, subdivision
plan or manifest references, output checksums, and artifact storage
location. Large POD5 artifacts must remain outside the package
repository.

## Status Semantics

Section and report-card statuses have intentionally narrow meanings:

- `pass`: evidence satisfies the configured QC threshold or requirement.
- `warn`: evidence is present but needs review or is incomplete/bounded.
- `fail`: evidence violates a configured QC threshold or requirement.
- `not_available`: the section is not part of this run or the
  integration has not landed yet.
- `not_checked`: the section is applicable but the check was not run.

Bamana validation failures that carry structured payloads remain QC
evidence. They should appear as `warn` or `fail` evidence instead of
being dropped as ordinary R exceptions. POD5 and library-preparation
evidence should follow the same principle once parser/runtime support is
complete.

## Synoptikon And Grammateus Boundaries

The payload is designed for Synoptikon ingestion and Grammateus report
construction, but it does not create a new identity, approval, or
signing system. Mneion remains the control-plane authority for governed
execution, storage context, evidence custody, authority keys, and audit
policy. Grammateus remains the authority for trusted technical report
manifests and lifecycle state.

## Validation

The floundeR test suite validates emitted JSON payloads against the
installed v1 JSON schema. It also checks for a canonical schema copy in
`../mnemosyne` or `../mnemosyne-docs` and compares it with the local
schema when present. No dedicated floundeR Synoptikon schema file exists
in those sibling repositories yet, so that sync check skips until the
control-plane copy is added.
