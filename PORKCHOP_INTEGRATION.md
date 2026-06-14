# Porkchop Integration Notes

This document records floundeR's curated Porkchop integration decisions.
Porkchop support must stay focused on nanopore library-preparation QC, review,
provenance, reporting, and synoptikon handoff. It must not turn floundeR into a
general trimming or preprocessing frontend.

## Audit Preflight

Checked on 2026-06-14:

- `../porkchop/AGENTS.md`
- `../porkchop/ROADMAP.md`
- `../porkchop/README.md`
- `../porkchop/docs/output/json.rst`
- `../porkchop/docs/output/schemas/`
- `../porkchop/docs/kits/provenance.rst`
- `../porkchop/docs/validation/index.rst`
- `../porkchop/docs/algorithms/scoring.rst`
- `../porkchop/docs/workflows/screening.rst`
- `../porkchop/docs/workflows/cdna.rst`
- `../porkchop/docs/limitations.rst`
- `../porkchop/src/lib.rs`
- `../porkchop/src/kit.rs`
- `../porkchop/src/motif_index.rs`
- `../porkchop/src/screen.rs`
- `../porkchop/src/cdna.rs`
- `../porkchop/src/trim.rs`

No Porkchop-specific canonical files were found in `../mnemosyne-docs` during
the floundeR governance review, so the local Porkchop agent instructions,
roadmap, README, Sphinx output/provenance/validation docs, public schemas, and
Rust source remain the governing boundary until canonical mnemosyne-docs
material exists.

Local Porkchop source state at audit time:

- crate: `porkchop`
- crate version: `0.3.234`
- git commit: `e60bdd2`
- MSRV: Rust `1.82`

floundeR is not linked to Porkchop yet. `flounder_rust_capabilities()` correctly
reports `porkchop = "not_linked"` until the curated bindings land.

## Public Contract Shape

Porkchop is a Rust toolkit for ONT kit inspection, adapter/barcode screening,
trimming, cDNA workflow support, and benchmarking. floundeR should consume it
as an in-process Rust library only where the output is QC evidence. The normal
R API must return data frames/lists and typed conditions rather than rendered
CLI text, terminal dashboards, HTML reports, or generic command passthroughs.

The public output contract most relevant to floundeR is the documented JSON
schema version `1.0.0` for screening and cDNA classification payloads. The
important semantic rule is that Porkchop's `score` and `normalized_score`
fields are heuristic evidence scores. They are not calibrated probabilities and
must not be described as posterior confidence, biological certainty, or
demultiplexing truth.

The kit registry is also a public contract for floundeR. Kit rows carry:

- kit id and description;
- chemistry and workflow;
- lifecycle status: `current`, `legacy`, `retired`, or `unknown`;
- support level: `full`, `best-effort`, `partial`, or `experimental`;
- introduced/retired years where known;
- kit and sequence provenance;
- source URLs;
- adapter, primer, flank, and barcode records.

Registry support and workflow validation are separate claims. A kit can be
present in the registry without being validated for trimming, demultiplexing,
or cDNA biological full-length behaviour. floundeR report outputs must carry
that distinction.

## First Curated QC Surface

The first floundeR Porkchop surface should promote only evidence-producing,
read-only contracts.

| floundeR target | Porkchop source | Purpose |
| --- | --- | --- |
| `library_kit_candidates()` | `list_supported_kits()`, `get_sequences_for_kit()`, kit metadata, and a future pure screening-ranker library function | Rank candidate library-preparation contexts from observed motif evidence while carrying lifecycle, support level, provenance, score terminology, and limitations. |
| `library_adapter_primer_evidence()` | `motif_index::cached_motif_index_for_kit()`, `KitMotifIndex`, `MotifFamily::Adapter`, `MotifFamily::Primer` | Report adapter and primer motif observations per read or aggregate sample without trimming reads. |
| `library_barcode_evidence()` | `KitMotifIndex`, `MotifFamily::Barcode`, `MotifFamily::Flank`, registry expected barcode counts | Report barcode/flank motif evidence, barcode ambiguity, expected barcode coverage, and barcode-family context. |
| `library_cdna_primer_evidence()` | `cdna::orientation_rules_for_kit()`, `cdna::detect_cdna_primer_pair()` | Report cDNA primer-pair evidence and read-class vocabulary for supported PCS/PCB kits. |
| `library_preparation_report_card()` | floundeR wrapper over the above evidence tables | Produce pass/warn/fail checks for unexpected library chemistry, adapter burden, barcode ambiguity, cDNA partial/unclassified burden, and unsupported or partial kit support levels. |

`library_kit_candidates()` needs the most care. Porkchop's current
`screen::run_screen()` is still a CLI-style workflow that owns input reading,
parallel execution, terminal UI, JSON writing, and optional HTML writing. Its
internal kit-ranking logic is useful, but it is not yet a clean floundeR
binding point because key helpers are private and output-oriented. Before
floundeR binds kit candidate screening, Porkchop should expose a pure library
adapter that accepts in-memory motif tallies or read batches and returns
structured candidate rows without writing files or starting the terminal UI.

`library_adapter_primer_evidence()` and `library_barcode_evidence()` can start
from the existing `motif_index` API. `KitMotifIndex` already preserves registry
order, motif family, source bucket, normalized motif text, reverse-complement
motifs, provenance through `SequenceRecord`, and ambiguity-aware match ranges.
The floundeR wrapper should convert those hits to stable R-native tables rather
than copying the CLI screen report shape wholesale.

`library_cdna_primer_evidence()` can start from the existing `cdna` API.
`detect_cdna_primer_pair()` takes a kit id and in-memory read sequence, returns
structured `CdnaPrimerPairDetection`, and distinguishes unsupported kit IDs
from supported kits with unclassified reads. floundeR should preserve that
distinction in R conditions or status columns.

## R Return-Shape Guidance

The first R wrappers should use stable, schema-versioned objects.

Kit candidate rows should include:

- `schema_version`;
- `kit_id`;
- `description`;
- `chemistry`;
- `workflow`;
- `kit_family`;
- `lifecycle_status`;
- `support_level`;
- `legacy`;
- `introduced_year`;
- `retired_year`;
- `score`;
- `normalized_score`;
- `score_kind = "heuristic_evidence_score"`;
- `matched_motifs`;
- `total_hits`;
- `provenance_source`;
- `provenance_appendix`;
- `provenance_notes`;
- `source_urls`;
- `validation_status`;
- `known_limitations`;

Motif evidence rows should include:

- `schema_version`;
- `read_id` or aggregate group id;
- `kit_id`;
- `motif_name`;
- `motif_kind`;
- `motif_family`;
- `motif_source`;
- `strand`;
- `start`;
- `end`;
- `match_semantics`;
- `edit_distance` when available;
- `provenance_source`;
- `support_level`;

cDNA evidence rows should include:

- `schema_version`;
- `kit_id`;
- `read_id`;
- `class`;
- `classified`;
- `full_length`;
- `five_prime_name`;
- `five_prime_start`;
- `five_prime_end`;
- `three_prime_name`;
- `three_prime_start`;
- `three_prime_end`;
- `primer_hit_count`;
- `workflow_support_note`;
- `known_limitations`;

Aggregate evidence and report-card helpers should preserve the same score and
support-level terminology. A bounded or sampled screen must identify the sample
fraction, sampled-read count, and evidence scope so it is not interpreted as a
full-file claim.

## Report-Card Checks

The initial report-card checks should be:

- unexpected library chemistry: warn/fail when the top supported kit candidate
  conflicts with supplied run metadata or expected kit id;
- adapter burden: warn/fail on adapter motif observations outside expected
  terminal contexts or above a configured aggregate burden;
- barcode ambiguity: warn/fail on tied barcode-family evidence, excessive
  flank-only evidence, or observed barcode counts inconsistent with the
  expected kit;
- cDNA partial/unclassified burden: warn/fail on high partial, primer-only, or
  unclassified rates when a cDNA workflow is expected;
- unsupported or partial support level: warn/fail when selected kit candidates
  have `partial`, `best-effort`, `experimental`, `unknown`, retired, or
  unsupported status.

These checks are QC and review signals. They should not perform trimming,
demultiplexing, read rescue, or biological full-length certification.

## Transformation Boundary

Porkchop contains trimming and cleaning workflows that are valuable in the
Porkchop project but out of scope for normal floundeR APIs. floundeR may report
motif evidence that would explain why trimming may be appropriate, but it
should not silently transform reads.

| Porkchop operation | floundeR posture | Boundary |
| --- | --- | --- |
| kit listing and description | read-only candidate | Acceptable as metadata/provenance for reports and synoptikon payloads. |
| screening/ranking | read-only candidate after pure library adapter | Acceptable as heuristic evidence only, never calibrated probability. |
| motif index matching | read-only candidate | Acceptable for adapter, primer, barcode, flank, and ambiguity evidence. |
| cDNA primer-pair detection | read-only candidate | Acceptable for cDNA evidence and class vocabulary; not full rescue/demultiplexing validation. |
| JSON schema/fixture validation | contract reference | Use schemas to shape R return values and tests; do not emit raw CLI JSON as the R API. |
| `clean` and trimming | out of normal floundeR API | Writes modified reads and belongs to preprocessing. Add only as explicit transformation workflow if a later roadmap slice requires it. |
| trimming coordinate decisions | explanatory evidence only for now | May inform QC findings, but not an implicit write-capable helper. |
| benchmarking | out of normal floundeR API | Engineering/performance evidence, not routine run QC. |
| terminal dashboard and HTML screen report | out of normal floundeR API | floundeR reports should go through Grammateus semantic elements. |
| whole-CLI passthrough | out of scope | floundeR must remain R-native and curated, not a generic Porkchop launcher. |

If a future write-capable helper is approved, it must require explicit output
paths, record input/output checksums, keep Porkchop/floundeR versions, preserve
the exact kit and threshold settings, and return a provenance manifest distinct
from read-only QC evidence.

## Report And Synoptikon Use

Porkchop-derived report sections should include:

- kit candidate table with score terminology and support-level caveats;
- adapter and primer evidence summaries;
- barcode/flank evidence summaries and ambiguity notes;
- cDNA primer-pair evidence and class distribution where cDNA workflows apply;
- registry provenance, lifecycle status, support level, and validation
  limitations;
- report-card findings for unexpected chemistry, adapter burden, barcode
  ambiguity, cDNA partial/unclassified burden, and unsupported support levels.

These outputs should feed future Grammateus report elements and the synoptikon
QC payload. Core library-preparation QC must remain available without private
Grammateus runtime assets.

## Required Upstream Library Work

No Porkchop source change is required for this audit. The next implementation
slice should either bind existing pure APIs directly or make small upstream
Porkchop changes to expose a clean screening adapter.

The likely upstream additions are:

- a public `ScreenReport`/`ScreenKitCandidate` data model in `porkchop::screen`;
- a pure function that accepts motif tallies or read batches and returns
  candidate rows, contexts, read counts, and sampling metadata;
- public helpers for canonical barcode naming and motif score weighting, if
  those remain part of the candidate score contract;
- optional serialization to the existing documented JSON schema without
  requiring the CLI `run_screen()` path.

Those additions must preserve Porkchop's own Sphinx docs, JSON schemas,
snapshot tests, score terminology, provenance requirements, and semantic
versioning expectations.
