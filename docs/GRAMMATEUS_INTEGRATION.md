# Grammateus Integration Notes

This document records floundeR’s Grammateus reporting integration
preflight and the boundaries that keep the public package open-source
while Grammateus remains private reporting infrastructure.

## Checkout Preflight

Checked on 2026-06-14 before adding any Rust dependency or R reporting
API:

- local checkout: `../grammateus`
- branch: `main`
- package version declared by Grammateus: `0.6.0`
- local commit: `0fba9ab3f23018133e4b52af47378e5a69519eff`
- remote: `origin` at `https://github.com/sagrudd/grammateus`
- remote comparison after `git fetch --prune`: `0` commits ahead and `0`
  commits behind `origin/main`
- worktree state: clean

This confirms the local private source checkout is current and clean
enough for the next Slice 13 steps: canonical contract verification and
public library API identification.

No Grammateus source code is vendored into floundeR by this preflight.
No floundeR build path, package check, or core QC API depends on private
Grammateus source or runtime assets at this stage.

## Distribution Boundary

The public floundeR package must continue to install, check, and provide
core QC functionality without Grammateus source. Development may use an
authorized private checkout to design and test in-process bindings, but
public distribution must use optional prebuilt Grammateus
runtime/reporting artifacts with explicit discovery, version validation,
manifests, and checksums.

The next integration steps must preserve the canonical Grammateus
contracts listed in `GOVERNANCE.md`, especially rendering elements,
trusted-report lifecycle, and trusted-report Synoptikon contracts.

## Canonical Contract Verification

Checked on 2026-06-14 before changing `../grammateus` or adding any
floundeR binding:

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

The documents are compatible with the floundeR reporting direction,
provided floundeR treats Grammateus as internal rendering/report
formalisation infrastructure rather than a user-facing product. The
binding and R API work must preserve these constraints:

- Do not expose `Grammateus` as customer-visible branding in reports or
  public workflow language. Customer-visible output may use Mnemosyne
  Biosciences branding supplied by the calling product/template.
- Do not make Grammateus responsible for nanopore QC decisions, report
  approval, user identity, role membership, authority-key rotation,
  storage, or workflow orchestration.
- Represent floundeR QC outputs as semantic tables, figures, plots,
  images, provenance, methods, limitations, and appendices with stable
  lower-snake-case identifiers, required captions, image alt text,
  deterministic rendering, and source/producers/version/timestamp
  provenance.
- Preserve the rendering-elements contract for controlled R/ggplot2
  execution: deterministic run directories, plot spec/data JSON,
  generated scripts, stdout/stderr, backend session metadata, artifact
  hashes, and SVG/PNG outputs.
- Reuse Grammateus report manifest, lifecycle, provenance-table, and
  verification semantics instead of creating a parallel trusted-report
  stack in floundeR.
- Reuse Synoptikon/Mneion managed users, domain roles, authority keys,
  governance-envelope signing, workbench policy, and GUI/API enforcement
  for trusted report lifecycle operations.
- Keep public floundeR install/check paths independent of private
  Grammateus source and private runtime downloads; private development
  bindings must be translated into the optional prebuilt runtime model
  before public release.

No conflict was found between the verified canonical documents and the
planned floundeR QC reporting interface. The next implementation slice
may identify the specific public Grammateus Rust APIs needed for report
elements, rendering, manifests, provenance, branding hooks, and
trusted-report lifecycle metadata.

## Public API Inventory For floundeR

Checked on 2026-06-14 against `../grammateus` commit
`0fba9ab3f23018133e4b52af47378e5a69519eff`.

The floundeR binding should stay narrow: expose report semantics needed
for nanopore QC, review, provenance, reporting, and Synoptikon handoff
without turning floundeR into a general Grammateus frontend. The current
Grammateus public Rust exports cover most of the required surface.

### Runtime And Version Metadata

Use these constants to record the private runtime/build used for
governed report rendering:

- `GRAMMATEUS_RELEASE`
- `ENGINE_VERSION`
- `TRUSTED_REPORT_MANIFEST_SCHEMA_VERSION`
- `REPORT_MANIFEST_GOVERNANCE_MESSAGE_TYPE`

floundeR manifests should record these alongside the floundeR, R,
`pod5-tools`, `bamana`, and `porkchop` versions used to generate a
report.

### Semantic Tables

Use Grammateus tables for run metadata, QC report cards, POD5 integrity
summaries, BAM summaries, Porkchop library-preparation evidence,
methods, limitations, and provenance appendices.

Required public exports:

- `ReportTable`
- `ReportTableColumn`
- `ReportTableValue`
- `ReportTableDataType`
- `ReportTableAlignment`
- `ReportTableLayoutHint`
- `ReportTableNote`
- `ReportElementProvenance`
- `ReportTable::validate()`
- `render_report_table_html()`
- `render_report_table_pdf()`

floundeR R helpers should convert data frames and lists into this model,
then validate before rendering. Table identifiers must keep the
Grammateus `table_` prefix and stable lower-snake-case form.

### Figures And Existing R Plot Artifacts

Use Grammateus figures when R has already produced an image artifact,
including `ggplot2` plots saved by floundeR as deterministic SVG/PNG
files.

Required public exports:

- `ReportFigure`
- `ReportFigureSource`
- `ReportImageSource`
- `ReportImageFormat`
- `ReportFigureLayoutHint`
- `ReportFigure::validate()`
- `render_report_figure_html()`
- `render_report_figure_pdf()`

This is the Rust-side target for
[`grammateus_figure_from_file()`](https://sagrudd.github.io/floundeR/reference/grammateus_figure_from_file.md)
and
[`grammateus_figure_from_ggplot()`](https://sagrudd.github.io/floundeR/reference/grammateus_figure_from_file.md).
floundeR must supply captions, alt text, methods notes where
appropriate, dimensions, checksums, and `ReportElementProvenance`.
Figure identifiers must keep the `figure_` prefix.

### Semantic Plots And Controlled R Execution

Use Grammateus plots when the report specification owns plot generation
and Grammateus calls a controlled R/ggplot2 backend.

Required public exports:

- `ReportPlot`
- `ReportPlotType`
- `ReportPlotData`
- `ReportPlotMappings`
- `ReportPlotAxes`
- `ReportPlotAxis`
- `ReportPlotAxisScale`
- `ReportBarValueSemantics`
- `ReportPlotOutput`
- `ReportPlotOutputFormat`
- `ReportPlot::validate()`
- `RPlotBackendOptions`
- `RPlotBackendExecution`
- `RPlotBackendRun`
- `RPlotBackendArtifact`
- `RPlotBackendError`
- `run_r_plot_backend()`

This is the Rust-side target for
[`grammateus_plot_spec()`](https://sagrudd.github.io/floundeR/reference/grammateus_plot_spec.md)
and
[`grammateus_render_plot()`](https://sagrudd.github.io/floundeR/reference/grammateus_plot_spec.md).
The first floundeR plot families should be mapped to the existing
`line`, `bar`, `stacked_bar`, and `scatter` plot families where
possible. Bar and stacked-bar QC plots must declare
`ReportBarValueSemantics` so counts, percentages, rates, and normalised
measurements are not confused.

### HTML/PDF Rendering

The current public rendering exports useful to floundeR are:

- `MarkdownReportOptions`
- `ReportMetadataTable`
- `ReportMetadataField`
- `OutputFormat`
- `render_markdown_report_html()`
- `render_markdown_report_pdf()`
- `render_output()`

These provide deterministic Mnemosyne-branded Markdown report rendering
and metadata title panels. They are suitable for the first transitional
[`qc_report()`](https://sagrudd.github.io/floundeR/reference/qc_report.md)
implementation if floundeR composes semantic elements into a canonical
Markdown payload with generated tables, figures, provenance, and
appendices.

Do not treat `render_rendering_elements_example_report()` or
`render_trusted_report_example_report()` as normal floundeR APIs. They
are useful as contract examples and test oracles only.

### Trusted Report Manifest, Provenance, And Lifecycle

Use Grammateus trusted-report primitives to avoid a parallel report
trust stack in floundeR.

Required public exports for manifests and artifact hashes:

- `ReportManifestArtifactHashes`
- `ReportManifestSigningAlgorithm`
- `ReportManifestSigningContext`
- `ReportSignedManifest`
- `ReportSignedManifest::new()`
- `validate_report_signed_manifest()`
- `canonical_serialize_report_signed_manifest()`
- `compute_report_signed_manifest_hash()`

Required public exports for provenance tables:

- `ReportTrustProvenanceRow`
- `ReportTrustProvenanceSignatureStatus`
- `build_report_trust_provenance_table()`

Required public exports for lifecycle metadata and validation:

- `ReportLifecycleState`
- `ReportLifecycleTransition`
- `ReportLifecycleRole`
- `ReportLifecycleActor`
- `ReportLifecyclePolicy`
- `ReportLifecycleTransitionRule`
- `ReportLifecycleEvent`
- `ReportLifecycleValidationError`
- `validate_report_lifecycle_transition()`

Required public exports for user attestations and governance envelopes:

- `ReportUserSignatureEvent`
- `ReportUserSignaturePurpose`
- `compute_report_user_signature_hash()`
- `validate_report_user_signature_event()`
- `validate_report_user_signature_chain()`
- `ReportAuthorityPublicKey`
- `ReportAuthorityKeyStatus`
- `ReportAuthorityKeyUsage`
- `ReportAuthorityKeyPublicFormat`
- `ReportGovernanceEnvelopeV1`
- `ReportGovernanceEnvelopeSigner`
- `ReportGovernanceEnvelopeVerification`
- `ReportGovernanceEnvelopeVerificationOutcome`
- `build_report_manifest_governance_envelope()`
- `canonical_serialize_report_governance_payload()`
- `compute_report_governance_envelope_hash()`
- `verify_report_manifest_governance_envelope()`

floundeR should construct, validate, and record report trust payloads,
but it must not own users, roles, authority-key rotation, final approval
decisions, or trusted-report orchestration. Synoptikon/Mneion remain
responsible for those governance actions.

### Branding And Templates

Use Grammateus-controlled branding and template hooks rather than
bespoke floundeR CSS or PDF styling.

Relevant public exports:

- `MarkdownReportOptions`
- `ReportMetadataTable`
- `ReportMetadataField`
- `Template`
- `TemplateDefinition`
- `TemplateRegistry`
- `TemplateGraph`

The current rendering functions already emit Mnemosyne Biosciences
branded HTML/PDF output for Markdown reports. Future floundeR bindings
should prefer a stable Grammateus theme/template profile when one is
exposed for full multi-element reports. Until then, floundeR should only
select approved Grammateus profiles such as the semantic plot
`theme`/`palette` fields and should not maintain separate report
branding assets.

### APIs Intentionally Out Of Scope

Do not bind these as normal floundeR APIs unless a later QC/reporting
slice identifies a direct need:

- Grammateus CLI/example report helpers beyond contract tests.
- Template authoring, qualification, approval, registry signing, and
  transparency-log administration.
- Generic document compilation, issuance, co-signing, bundle export, and
  quotation/document DSL workflows.
- Any API that would make floundeR responsible for Synoptikon/Mneion
  identity, authority keys, governance policy enforcement, or report
  lifecycle orchestration.

### API Gaps To Track

The public exports are sufficient for the next R-to-Grammateus figure
handoff and controlled R plot-generation wrappers. Before the final
[`qc_report()`](https://sagrudd.github.io/floundeR/reference/qc_report.md)
surface is considered complete, floundeR should either use or request
these stable Grammateus library APIs:

- a multi-element report composer that accepts ordered `ReportTable`,
  `ReportFigure`, `ReportPlot`, methods, limitations, appendix, and
  provenance elements and emits HTML/PDF plus manifest inputs without
  first reducing the report to Markdown;
- a stable public theme/profile selector for Mnemosyne-branded QC
  reports;
- a public runtime bundle metadata API for prebuilt private Grammateus
  assets, including version, platform, ABI, manifest, checksum, and
  compatibility claims;
- JSON/serde constructors or validation helpers that make the R
  extension boundary predictable for R lists and data frames.
