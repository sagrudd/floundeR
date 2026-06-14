# Grammateus Integration Notes

This document records floundeR's Grammateus reporting integration preflight and
the boundaries that keep the public package open-source while Grammateus remains
private reporting infrastructure.

## Checkout Preflight

Checked on 2026-06-14 before adding any Rust dependency or R reporting API:

- local checkout: `../grammateus`
- branch: `main`
- package version declared by Grammateus: `0.6.0`
- local commit: `0fba9ab3f23018133e4b52af47378e5a69519eff`
- remote: `origin` at `https://github.com/sagrudd/grammateus`
- remote comparison after `git fetch --prune`: `0` commits ahead and `0`
  commits behind `origin/main`
- worktree state: clean

This confirms the local private source checkout is current and clean enough for
the next Slice 13 steps: canonical contract verification and public library API
identification.

No Grammateus source code is vendored into floundeR by this preflight. No
floundeR build path, package check, or core QC API depends on private
Grammateus source or runtime assets at this stage.

## Distribution Boundary

The public floundeR package must continue to install, check, and provide core
QC functionality without Grammateus source. Development may use an authorized
private checkout to design and test in-process bindings, but public
distribution must use optional prebuilt Grammateus runtime/reporting artifacts
with explicit discovery, version validation, manifests, and checksums.

The next integration steps must preserve the canonical Grammateus contracts
listed in `GOVERNANCE.md`, especially rendering elements, trusted-report
lifecycle, and trusted-report Synoptikon contracts.

## Canonical Contract Verification

Checked on 2026-06-14 before changing `../grammateus` or adding any floundeR
binding:

- `../mnemosyne-docs/charters/products/grammateus/GRAMMATEUS_CHARTER.md`
- `../mnemosyne-docs/requirements/srs/products/grammateus/GRAMMATEUS_SRS.md`
- `../mnemosyne-docs/architecture/ards/products/grammateus/ARD-001-grammateus-architecture-review.md`
- `../mnemosyne-docs/architecture/ards/products/grammateus/ARD-002-quotation-document-dsl-architecture-review.md`
- `../mnemosyne-docs/architecture/ards/products/grammateus/SOURCE_MAP.md`
- `../mnemosyne-docs/products/grammateus/GRAMMATEUS_RENDERING_ELEMENTS_CONTRACTS.md`
- `../mnemosyne-docs/products/grammateus/GRAMMATEUS_TRUSTED_REPORT_LIFECYCLE_CONTRACTS.md`
- `../mnemosyne-docs/products/grammateus/GRAMMATEUS_TRUSTED_REPORT_SYNOPTIKON_CONTRACTS.md`
- `../mnemosyne-docs/products/grammateus/README.md`
- `../mnemosyne-docs/products/grammateus/TAXONOMY_MAP.md`

The documents are compatible with the floundeR reporting direction, provided
floundeR treats Grammateus as internal rendering/report formalisation
infrastructure rather than a user-facing product. The binding and R API work
must preserve these constraints:

- Do not expose `Grammateus` as customer-visible branding in reports or public
  workflow language. Customer-visible output may use Mnemosyne Biosciences
  branding supplied by the calling product/template.
- Do not make Grammateus responsible for nanopore QC decisions, report approval,
  user identity, role membership, authority-key rotation, storage, or workflow
  orchestration.
- Represent floundeR QC outputs as semantic tables, figures, plots, images,
  provenance, methods, limitations, and appendices with stable lower-snake-case
  identifiers, required captions, image alt text, deterministic rendering, and
  source/producers/version/timestamp provenance.
- Preserve the rendering-elements contract for controlled R/ggplot2 execution:
  deterministic run directories, plot spec/data JSON, generated scripts,
  stdout/stderr, backend session metadata, artifact hashes, and SVG/PNG outputs.
- Reuse Grammateus report manifest, lifecycle, provenance-table, and
  verification semantics instead of creating a parallel trusted-report stack in
  floundeR.
- Reuse Synoptikon/Mneion managed users, domain roles, authority keys,
  governance-envelope signing, workbench policy, and GUI/API enforcement for
  trusted report lifecycle operations.
- Keep public floundeR install/check paths independent of private Grammateus
  source and private runtime downloads; private development bindings must be
  translated into the optional prebuilt runtime model before public release.

No conflict was found between the verified canonical documents and the planned
floundeR QC reporting interface. The next implementation slice may identify
the specific public Grammateus Rust APIs needed for report elements, rendering,
manifests, provenance, branding hooks, and trusted-report lifecycle metadata.
