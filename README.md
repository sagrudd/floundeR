# floundeR

`floundeR` is an R-native toolbox for Oxford Nanopore sequence-data QC and
review. It reads the run products analysts already use in R, returns stable
tidy QC contracts, and uses in-process Rust where the heavy lifting belongs:
POD5 raw-data inspection through `pod5-tools`, BAM/BGZF/FASTQ evidence through
`bamana`, and library-preparation evidence through `porkchop`.

<!-- badges: start -->
[![Code size](https://img.shields.io/github/languages/code-size/sagrudd/floundeR)](https://img.shields.io/github/languages/code-size/sagrudd/floundeR)
[![Releases](https://img.shields.io/github/downloads/sagrudd/floundeR/total)](https://img.shields.io/github/downloads/sagrudd/floundeR/total)
<!-- badges: end -->

The reboot goal is not to mirror every command from the Rust projects. The goal
is to make `floundeR` the best R interface for nanopore run health, raw-data
integrity, alignment-level QC, library-preparation evidence, provenance,
technical reporting, and Synoptikon handoff.

FAST5 support has been retired. Current raw-signal QC is POD5-oriented.

## What floundeR Does Now

- Summarises MinKNOW, Guppy, and Dorado sequencing-summary files into
  schema-versioned tidy QC tables.
- Builds run-level report cards for yield, read counts, pass fraction, Q-score,
  read length, channel density, and barcode balance.
- Discovers, verifies, manifests, compares, and plans demonstration subdivision
  for POD5 folders without shelling out to a CLI.
- Uses curated Bamana-backed BAM checks for summaries, verification,
  validation, EOF/index/map/sort/tag evidence, and BAM report cards.
- Uses curated Porkchop-backed evidence for ONT kit candidates, adapter and
  primer motifs, barcode evidence, cDNA primer evidence, and
  library-preparation report cards.
- Exports versioned Synoptikon QC payloads.
- Builds Grammateus-shaped report contracts, governed figures, plot specs, and
  manifests while keeping the private Grammateus runtime optional.

## Installation

Install source-build requirements first:

- R with the declared CRAN/Bioconductor package dependencies.
- Cargo and `rustc` matching the `SystemRequirements` floor in `DESCRIPTION`.
- Platform build tools for compiling an R package with native code.

Then install from GitHub:

```r
install.packages("devtools")
devtools::install_github("sagrudd/floundeR")
```

Detailed source-install notes for macOS, Linux, Windows, Docker, and CI are in
[`DEVELOPMENT.md`](DEVELOPMENT.md#source-install-requirements). Public GitHub
installation paths and optional private Grammateus runtime setup are documented
in [`GITHUB_INSTALLATION.md`](GITHUB_INSTALLATION.md).

Core QC does not require private Grammateus source or runtime assets. Governed
HTML/PDF rendering can use an optional prebuilt Grammateus runtime when
authorized and configured through `GRAMMATEUS_HOME`, package options, or the
runtime cache helpers.

## Run Health In One Minute

This is the smallest useful health check for a nanopore run. It uses the
packaged sequencing-summary fixture here, but the same code accepts a real
MinKNOW, Guppy, or Dorado sequencing-summary path.

```r
library(floundeR)

summary_file <- flnDr("sequencing_summary.txt.bz2")

run_summary <- qc_run_summary(summary_file)
barcode_balance <- qc_barcode_composition(summary_file)
run_card <- qc_report_card(
  run_summary,
  barcode_composition = barcode_balance
)

run_summary[, c(
  "read_count",
  "total_bases",
  "pass_fraction",
  "mean_qscore",
  "n50_read_length",
  "channel_count"
)]

table(run_card$status)

payload <- as_synoptikon_qc(
  run_summary = run_summary,
  barcode = barcode_balance,
  report_cards = list(run = run_card)
)
payload$schema_version
```

That gives an analyst a one-row run summary, pass/warn/fail report-card counts,
and a versioned payload shape ready for downstream handoff.

## Contemporary QC Workflow

Start with the sequencing-summary file. `floundeR` accepts a path, a
`SequencingSummary` object, or a normalised data frame/tibble for the core QC
helpers.

```r
library(floundeR)

summary_file <- flnDr("sequencing_summary.txt.bz2")

run_summary <- qc_run_summary(summary_file)
yield_over_time <- qc_yield_over_time(summary_file, resolution_minutes = 15)
length_distribution <- qc_read_length_distribution(summary_file)
quality_distribution <- qc_quality_distribution(summary_file)
channel_density <- qc_channel_density(summary_file)
barcode_balance <- qc_barcode_composition(summary_file)

run_card <- qc_report_card(
  run_summary,
  barcode_composition = barcode_balance
)
```

The outputs are ordinary tibbles/data frames with explicit schema versions, so
they can be inspected interactively, plotted, placed in reports, or exported to
downstream systems.

## POD5 Raw-Data Evidence

POD5 helpers call the embedded Rust extension and return R-native objects. They
are intended for raw-data inventory, integrity checks, provenance, and report
evidence.

```r
run_dir <- "/path/to/nanopore/run"

pod5_folders <- pod5_find(run_dir)
pod5_inventory <- pod5_folder_info(run_dir)
pod5_files <- pod5_manifest(run_dir)
```

For a single local POD5 file:

```r
pod5_verify("/path/to/read_batch.pod5")
pod5_file_info("/path/to/read_batch.pod5")
```

`floundeR` also records one canonical opt-in real-data POD5 source from ONT
open data. New examples, vignettes, integration checks, and report demos should
start from these helpers rather than choosing arbitrary objects from the public
bucket:

```r
ont_zymo_pod5_dataset()
ont_zymo_pod5_example_objects()
```

The selected routine example is
`PAU85136_pass_279c9095_68316534_8289.pod5` from
`s3://ont-open-data/zymo_fecal_2025.05/raw/PAU85136/pod5/`. The selected
fail-state example is `PAU85136_fail_279c9095_68316534_0.pod5`. These files are
large enough that downloads are always explicit and opt-in:

```r
if (identical(tolower(Sys.getenv("FLOUNDER_RUN_NETWORK_TESTS")), "true")) {
  object <- ont_zymo_pod5_example_objects(role = "pass")
  ont_open_data_fetch(
    key = object$key,
    cache_dir = tools::R_user_dir("floundeR", "cache")
  )
}
```

Downloaded POD5 files and large derived POD5 artifacts must not be committed to
this repository. Dataset policy and provenance requirements are documented in
[`DATASETS.md`](DATASETS.md).

## BAM And Library-Preparation QC

Alignment evidence is handled through curated Bamana-backed functions:

```r
bam <- bam_summary("/path/to/alignments.bam")
bam_card <- bam_qc_report_card(
  summary = bam,
  eof = bam_check_eof("/path/to/alignments.bam"),
  index = bam_check_index("/path/to/alignments.bam"),
  sort = bam_check_sort("/path/to/alignments.bam")
)
```

Library-preparation evidence is handled through curated Porkchop-backed
functions. Scores are heuristic evidence, not calibrated probabilities.

```r
reads <- c(
  "AATGTACTTCGTTCAGTTACGTATTGCT",
  "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT"
)

kit_candidates <- library_kit_candidates(reads)
adapter_evidence <- library_adapter_primer_evidence(reads, kit_id = "LSK114")

library_card <- library_preparation_report_card(
  kit_candidates = kit_candidates,
  adapter_primer = adapter_evidence,
  expected_kit_id = "LSK114"
)
```

`floundeR` deliberately exposes only the subset of Bamana and Porkchop behavior
that improves nanopore QC, review, provenance, reporting, and Synoptikon
handoff.

## Synoptikon Handoff

Use the Synoptikon payload helpers when a QC result needs to move into the
Mnemosyne/Synoptikon lifecycle.

```r
payload <- as_synoptikon_qc(
  run_summary = run_summary,
  barcode = barcode_balance,
  pod5 = pod5_files,
  bam = bam,
  library_preparation = kit_candidates,
  report_cards = list(
    run = run_card,
    bam = bam_card,
    library_preparation = library_card
  )
)

write_synoptikon_qc(
  tempfile(fileext = ".json"),
  run_summary = run_summary,
  report_cards = list(run = run_card)
)
```

The payload schema is versioned and intended for downstream validation.

## Grammateus Reporting

Production QC reports are moving away from RMarkdown-first workflows and toward
Grammateus semantic report contracts. Public `floundeR` builds can prepare
report contracts, figures, plot specs, manifests, and provenance without the
private Grammateus runtime.

R can hand an existing plot artifact into a governed report:

```r
plot_spec <- qc_plot_yield_over_time(yield_over_time)
plot_run <- grammateus_render_plot(
  plot_spec,
  execution = "local_rscript",
  run_root = tempdir()
)
figures <- lapply(plot_run$artifacts, `[[`, "figure")

theme <- grammateus_mnemosyne_theme()
elements <- grammateus_apply_theme(
  grammateus_qc_report_elements(
    qc_summary = run_summary,
    yield_over_time = yield_over_time,
    quality_distribution = quality_distribution,
    barcode_balance = barcode_balance,
    report_card_findings = run_card
  ),
  theme
)

report <- qc_report(
  elements = elements,
  figures = figures,
  output_dir = tempdir(),
  output = c("html", "pdf")
)
```

When the private runtime is absent, `qc_report()` still writes the report
contract and manifest. Runtime-backed HTML/PDF rendering is enabled only after
an authorized prebuilt Grammateus runtime has been installed and validated:

```r
grammateus_runtime_available()
grammateus_runtime_validate()
```

The detailed interface is documented in
[`REPORTING_INTERFACE.md`](REPORTING_INTERFACE.md), and the private-runtime
distribution boundary is documented in [`DISTRIBUTION.md`](DISTRIBUTION.md).

## Development And Governance

- GitHub Actions are currently disabled during reboot churn; templates live in
  `.github/disabled-workflows/`.
- Local dependency and container build guidance is in
  [`DEVELOPMENT.md`](DEVELOPMENT.md) and [`DEPENDENCIES.md`](DEPENDENCIES.md).
- Cross-repository boundaries are tracked in [`GOVERNANCE.md`](GOVERNANCE.md).
- Future Bioconductor submission posture is tracked in
  [`BIOCONDUCTOR_POLICY.md`](BIOCONDUCTOR_POLICY.md).

The package documentation site is available at
[https://sagrudd.github.io/floundeR/](https://sagrudd.github.io/floundeR/).
