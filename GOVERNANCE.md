# Governance And Cross-Repository Boundaries

This document records the first-pass governance boundaries for the floundeR
revival. It is intentionally practical: agents should inspect these files
before making cross-repository changes.

## floundeR

No canonical `floundeR` product charter was found in `../mnemosyne-docs` during
the 2026-06-13 baseline review. Until one is authored, this repository's
governing local files are:

- `AGENTS.md`
- `ROADMAP.md`
- `TODO.md`
- `DATASETS.md`
- `DISTRIBUTION.md`
- `GITHUB_INSTALLATION.md`
- `REPORTING_INTERFACE.md`
- `GRAMMATEUS_RUNTIME_INTERFACE.md`
- `GRAMMATEUS_RELEASE_ASSETS.md`
- `LEGACY_REPORTING.md`

## pod5-tools

No canonical `pod5-tools` product charter was found in `../mnemosyne-docs`
during the 2026-06-13 baseline review. Until one is authored, inspect:

- `../pod5-tools/AGENTS.md`
- `../pod5-tools/README.md`
- `../pod5-tools/roadmap.md`
- `../pod5-tools/todo.md`

Only the subset of `pod5-tools` needed for nanopore QC, review, provenance,
demonstration datasets, reporting, and synoptikon handoff should be exposed
through `floundeR`.

## Bamana

No canonical Bamana product charter was found in `../mnemosyne-docs` during the
2026-06-13 baseline review. Until one is authored, Bamana integration is
governed by:

- `../bamana/CHARTER.md`
- `../bamana/docs/project-charter.md`
- `../bamana/docs/json-output.md`
- Bamana's public schemas/examples where present

Only BAM/BGZF/FASTQ behavior that directly supports nanopore QC, review,
provenance, reporting, and synoptikon handoff should be surfaced in `floundeR`.

## Porkchop

No canonical Porkchop product charter was found in `../mnemosyne-docs` during
the 2026-06-13 baseline review. Until one is authored, Porkchop integration is
governed by:

- `../porkchop/AGENTS.md`
- `../porkchop/ROADMAP.md`
- `../porkchop/README.md`
- `../porkchop/docs/output/json.rst`
- `../porkchop/docs/kits/provenance.rst`
- `../porkchop/docs/validation/index.rst`

Only adapter, primer, barcode, kit-registry, cDNA, and library-preparation
evidence that directly supports nanopore QC, review, provenance, reporting, and
synoptikon handoff should be surfaced in `floundeR`.

## Grammateus

Grammateus has canonical documentation in `../mnemosyne-docs`. Inspect at
least:

- `../mnemosyne-docs/charters/products/grammateus/GRAMMATEUS_CHARTER.md`
- `../mnemosyne-docs/requirements/srs/products/grammateus/GRAMMATEUS_SRS.md`
- `../mnemosyne-docs/products/grammateus/GRAMMATEUS_RENDERING_ELEMENTS_CONTRACTS.md`
- `../mnemosyne-docs/products/grammateus/GRAMMATEUS_TRUSTED_REPORT_LIFECYCLE_CONTRACTS.md`
- `../mnemosyne-docs/products/grammateus/GRAMMATEUS_TRUSTED_REPORT_SYNOPTIKON_CONTRACTS.md`

Grammateus should provide report element, provenance, rendering, and trusted
report semantics. `floundeR` should not create a parallel report trust stack.
Grammateus source code is private. `floundeR` must therefore interact with
Grammateus through optional, prebuilt, verified runtime/reporting artifacts or
through authorized private development checkouts, not by making private source a
mandatory dependency of the open-source package.

## Synoptikon And Mneion

Synoptikon-facing work is expected to land through `../mnemosyne`. The relevant
canonical Mneion control-plane documents in `../mnemosyne-docs` include:

- `../mnemosyne-docs/charters/products/mneion/MNEION_CHARTER.md`
- `../mnemosyne-docs/requirements/srs/products/mneion/MNEION_SRS.md`
- `../mnemosyne-docs/architecture/ards/products/mneion/ARD-001-mneion-architecture-review.md`

Any QC payload, report lifecycle, signing, or governance workflow must preserve
those boundaries and reuse existing Synoptikon/Mneion identity, role,
authority-key, audit, and governance-signing patterns.

## Generated Documentation Policy

During normal development, source documentation should be edited first and
generated artifacts should be refreshed intentionally. The repository may keep
generated `man/` files because they are part of normal R package development.
Generated `docs/`/pkgdown output should be treated as release-oriented output:
update it when public documentation is intentionally refreshed, but do not let
stale generated pages drive source behavior.

Root-level project-control files such as `AGENTS.md`, `ROADMAP.md`, `TODO.md`,
`DATASETS.md`, `DISTRIBUTION.md`, `GITHUB_INSTALLATION.md`,
`REPORTING_INTERFACE.md`,
`GRAMMATEUS_RUNTIME_INTERFACE.md`, `GRAMMATEUS_RELEASE_ASSETS.md`,
`LEGACY_REPORTING.md`, and this file are repository governance files. They are
excluded from R source package builds through `.Rbuildignore`.
