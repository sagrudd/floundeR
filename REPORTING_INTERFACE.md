# Grammateus Reporting Interface

`floundeR` reports should be built through Grammateus, not as RMarkdown-first
documents. R remains central for statistics and plotting, but Grammateus owns
the report contract: semantic elements, rendering, provenance, manifests,
branding, and trusted-report lifecycle compatibility.

`floundeR` itself remains open-source. Grammateus remains private, so the public
package must treat Grammateus as an optional prebuilt runtime/reporting bundle.
Core QC APIs must remain usable when Grammateus is absent; only governed report
rendering should require the private runtime.

## Goals

- Let R users drop existing `ggplot2` plots into governed QC reports.
- Let Grammateus call controlled R plot generation when the report
  specification owns the build.
- Preserve Mnemosyne Biosciences branding and standard technical-report
  structure.
- Keep PDF and HTML rendering deterministic and provenance-rich.
- Avoid a second report trust stack in `floundeR`.

## Direction 1: R Calls Grammateus

This is the primary interactive R workflow.

Expected R-facing API shape:

```r
qc <- qc_run_summary(...)

plot <- qc_plot_yield(qc)

figure <- grammateus_figure_from_ggplot(
  plot = plot,
  figure_id = "figure_yield_over_time",
  caption = "Cumulative sequencing yield over time.",
  alt_text = "Line chart showing cumulative bases over elapsed sequencing time.",
  methods_note = "Generated from sequencing summary records.",
  output_dir = cache_dir
)

report <- qc_report(
  elements = grammateus_qc_report_elements(qc_summary = qc),
  figures = list(figure),
  output_dir = cache_dir,
  output = c("html", "pdf")
)
```

The R wrapper should:

- save the plot to deterministic SVG/PNG artifacts;
- compute artifact checksums;
- collect source-data provenance;
- construct a Grammateus `ReportFigure` or `ReportPlot` equivalent;
- call the embedded Grammateus Rust renderer through R with
  `grammateus_render_element()`, `grammateus_render_figure_html()`, or
  `grammateus_render_figure_pdf()`;
- return output paths plus a manifest/provenance payload.

Implemented floundeR integration:

`qc_report()` assembles prepared semantic elements and governed figures into a
`flounder.qc_report_contract.v1` JSON file and a
`flounder.qc_report_manifest.v1` JSON file. Public builds always produce those
artifacts without private Grammateus assets. Requested HTML/PDF outputs are
recorded in the manifest as `runtime_unavailable`, `not_requested`, or
`pending_runtime_binding` unless an authorized private runtime provides the
rendering path.

This path is best when an analyst is already working in R and wants to add
custom exploratory or domain-specific plots to a governed report.

## Direction 2: Grammateus Calls R

This is the governed report-build workflow.

Grammateus already exposes backend-neutral `ReportPlot` specifications and a
controlled R/ggplot2 backend that can run local `Rscript` or a Docker-compatible
container. It records deterministic run directories, plot specs, data JSON,
generated scripts, stdout/stderr, `sessionInfo()`, artifact hashes, and rendered
SVG/PNG files.

Implemented floundeR integration:

```r
plot_spec <- grammateus_plot_spec(
  plot_id = "plot_quality_distribution",
  plot_type = "bar",
  data = quality_bins,
  mappings = list(x = "qscore_bin", y = "read_count"),
  axes = list(
    x = list(label = "Mean Q-score bin"),
    y = list(label = "Reads")
  ),
  caption = "Distribution of mean read quality.",
  value_semantics = "counts"
)

run <- grammateus_render_plot(
  plot_spec,
  execution = "local_rscript",
  run_root = cache_dir
)
```

`grammateus_plot_spec()` returns a `flounder_grammateus_plot_spec` object that
matches the Grammateus `ReportPlot` contract closely enough for later Rust
runtime binding: stable plot identity, semantic plot family, inline tidy data or
data reference, mappings, axes, output settings, optional bar-value semantics,
and provenance. `grammateus_render_plot()` writes a deterministic run directory
containing `plot_spec.json`, `plot_data.json`, `render_plot.R`,
`backend_session.txt`, stdout/stderr capture, and requested PNG/SVG artifacts.
The generated artifact metadata includes SHA-256 hashes and
`flounder_grammateus_figure` wrappers so the outputs can be handed into the
reporting path consistently.

The first render-binding surface is now present in floundeR:
`grammateus_render_element()`, `grammateus_render_figure_html()`, and
`grammateus_render_figure_pdf()` call registered Rust functions inside the R
extension. Public builds do not link the private Grammateus renderer and return
a typed `flounder_grammateus_runtime_unavailable` condition. Authorized builds
can replace that Rust response path with private Grammateus library/runtime
calls without changing the R report-element API.

This path is best when a report template or synoptikon workflow owns the report
definition and wants reproducible plot generation from canonical tidy data.

## Required R Helpers

`floundeR` should provide a thin R interface over Grammateus rather than
requiring users to hand-author Rust-shaped JSON.

R helpers:

- `grammateus_plot_spec()`: build and validate a backend-neutral plot spec from
  tidy R data.
- `grammateus_figure_from_file()`: wrap a PNG/SVG file as a governed figure.
- `grammateus_figure_from_ggplot()`: save a `ggplot2` object and wrap the
  artifact as a governed figure.
- `grammateus_render_plot()`: ask Grammateus to run its controlled R backend for
  a semantic plot spec.
- `grammateus_render_element()`: render a prepared semantic element through the
  registered Rust report-rendering boundary.
- `grammateus_render_figure_html()` and `grammateus_render_figure_pdf()`:
  figure-specific HTML/PDF render bindings for report assembly.
- `grammateus_report_element()`: low-level conversion for advanced users.
- `qc_report()`: high-level floundeR report builder that emits HTML, PDF, and
  manifest/provenance outputs.

These helpers should be R-native and ergonomic, but the generated report
elements must remain Grammateus semantic elements.

## Required QC Plot Families

The first floundeR report interface should cover:

- `plot_yield_over_time`: cumulative bases or reads over elapsed sequencing
  time.
- `plot_quality_distribution`: read quality distribution.
- `plot_read_length_distribution`: read length distribution.
- `plot_flowcell_density`: channel density or spatial yield.
- `plot_barcode_balance`: barcode read/yield balance.
- `plot_pod5_integrity`: pass/warn/fail POD5 integrity summaries.
- `plot_bam_mapping_summary`: mapped/unmapped/primary/supplementary/secondary
  BAM categories.
- `plot_bam_mapq_distribution`: mapping quality distribution.
- `plot_bam_flag_summary`: duplicate, QC-fail, paired, read1/read2, and strand
  categories.

Specialized QC plots can still be supplied as `ggplot2` objects through
`grammateus_figure_from_ggplot()`.

## Element Rules

Every report element produced by floundeR must have:

- a stable lower-snake-case identifier with a type prefix such as
  `plot_yield_over_time` or `figure_flowcell_density`;
- a caption;
- alt text for image-like HTML output;
- source-data provenance;
- producer name and version;
- production timestamp;
- artifact checksum where a file is rendered;
- deterministic dimensions and output format.

Plot data should be tidy. When data are too large to inline, report elements
should reference a canonical data artifact with a checksum and row count.

## Branding And Technical Structure

Reports must use Mnemosyne Biosciences branding through Grammateus templates or
themes. `floundeR` should not maintain one-off CSS or parallel PDF styling.

Implemented floundeR integration:

```r
theme <- grammateus_mnemosyne_theme()
report <- grammateus_apply_theme(grammateus_qc_report_elements(...), theme)
```

`grammateus_mnemosyne_theme()` records the Mnemosyne Biosciences brand, approved
Grammateus template reference, approved theme reference, palette/profile,
status colours, style policy, and provenance. The helper is intentionally
runtime-free: it does not embed private Grammateus assets, inline CSS, or
RMarkdown styling. `grammateus_apply_theme()` attaches that descriptor to
prepared semantic elements so later `qc_report()` and authorized Grammateus
runtime builds can resolve the private templates/themes through the governed
runtime manifest.

The standard technical QC report structure should be:

1. Title page and run identity.
2. Executive QC report card.
3. Input data provenance.
4. Sequencing/run summary.
5. POD5 raw-data integrity and metadata findings where available.
6. BAM/alignment QC findings where available.
7. Library-preparation findings where available, including Porkchop-derived
   adapter, primer, barcode, kit, and cDNA evidence.
8. Barcode/sample findings where available.
9. Methods and software versions.
10. Limitations.
11. Appendices and provenance table.

## Enterprise Requirements

- The interface must work in local analyst mode and governed synoptikon mode.
- The open-source package must not require Grammateus source code to install or
  check.
- Grammateus runtime discovery, version checks, and manifest/signature
  validation must be explicit and actionable.
- Governed mode should prefer containerized R execution through Grammateus and
  Synoptikon/Mneion execution controls.
- R package versions, generated scripts, input data, output artifacts, stdout,
  stderr, and session metadata must be recorded.
- Report rendering should be reproducible enough for review, approval, and
  later verification.
- Failures should be structured R conditions, not unparsed process logs.

## Out Of Scope

- `floundeR` should not become a general report-rendering framework.
- `floundeR` should not bypass Grammateus manifest, provenance, or trusted
  report contracts.
- `floundeR` should not invent an independent identity, approval, or signing
  model.
- `floundeR` should not present Porkchop heuristic screening scores as
  calibrated probabilities.
