# floundeR 0.0.24

* Reconciled the current import surface by declaring narrow namespace imports
  for active FASTA, FASTQ, GenBank, plotting, and ONT open-data dependencies.
* Removed stale RMarkdown-era packages from the development dependency
  bootstrap list now that legacy vignettes are excluded from package checks.
* Documented the remaining legacy `%>%` and `BamFile` exports so release-style
  package checks can complete without warnings.

# floundeR 0.0.23

* Added a release-style package check helper that builds and checks a source
  tarball outside the repository, avoiding diagnostics from development-only
  files and local check output.
* Corrected source-build ignore rules so repository governance files,
  development scripts, `.github`, and `.dockerignore` are excluded from package
  tarballs.

# floundeR 0.0.22

* Moved legacy RMarkdown vignettes to repo-only historical reference material
  while executable report workflows move to the Grammateus migration.

# floundeR 0.0.21

* Corrected POD5 manifest fixture checksum values so they remain character
  SHA-256-like hex strings under package checks and schema tests.

# floundeR 0.0.20

* Added a small installed BLAST fixture and restored `Blast$count()` behavior
  so the BLAST example and test no longer depend on a missing large UniRef
  output file.

# floundeR 0.0.19

* Added conservative sequencing-summary column alias handling for current
  Dorado and MinKNOW/POD5-era summary shapes.
* Allowed reduced Dorado summary files without `passes_filtering` to import as
  partial QC inputs with an explicit warning and `NA` logical values.

# floundeR 0.0.18

* Replaced generated `eval(parse())` sequencing-summary column specifications
  with structured `readr::cols_only()` collector construction.
* Added regression coverage for legacy Guppy-like and Dorado-like sequencing
  summary fixtures.

# floundeR 0.0.17

* Added a Docker-backed development build path and a dependency audit script for
  assessing declared CRAN/Bioconductor package versions and support channels in
  a clean container.
* Recorded the current R 4.6.0/Bioconductor 3.23 dependency snapshot and the
  package-quality failures exposed once dependencies are available.

# floundeR 0.0.16

* Added fixture-backed schema tests for the initial sequencing-summary,
  barcode-summary, and POD5-manifest QC table contracts.

# floundeR 0.0.15

* Added a test helper and development guidance so optional compiled Rust-backed
  tests skip cleanly until in-process bindings are available.

# floundeR 0.0.14

* Added a test helper and development guidance that keep network-dependent ONT
  open-data/S3 tests opt-in through `FLOUNDER_RUN_NETWORK_TESTS=true` and
  skipped by default on CRAN-like checks.

# floundeR 0.0.13

* Added small offline test fixtures for sequencing summaries, barcoding
  summaries, FASTQ, FASTA, and POD5 metadata manifests without committing
  binary POD5 files.

# floundeR 0.0.12

* Added a dry-run-first R dependency bootstrap script and development notes for
  local checks after FAST5 retirement.
* Configured the package for `testthat` edition 3.

# floundeR 0.0.11

* Retired the active FAST5 API surface by removing the `Fast5` export,
  implementation, fixtures, tests, vignette, generated reference page, and
  `rhdf5` dependency. Raw-signal QC now points toward POD5-era integration.

# floundeR 0.0.10

* Added a repository-level governance check that enforces the open-source
  boundary for floundeR, pod5-tools, bamana, and porkchop while keeping
  Grammateus private, optional, and distributed through prebuilt runtime assets.

# floundeR 0.0.9

* Added Porkchop to the curated open-source Rust engine plan for adapter,
  primer, barcode, kit-registry, and cDNA/library-preparation QC evidence.

# floundeR 0.0.8

* Documented the distribution boundary: `floundeR`, `pod5-tools`, and `bamana`
  stay open-source, while Grammateus remains private and is consumed through
  optional prebuilt runtime/reporting artifacts.

# floundeR 0.0.7

* Added an explicit Grammateus reporting interface plan for dropping R plots
  into governed HTML/PDF reports and for controlled R plot generation from
  Grammateus report specifications.

# floundeR 0.0.6

* Started the reboot planning baseline for contemporary nanopore QC and review.
* Documented the FAST5 retirement direction, selective Rust-in-R integration
  with `pod5-tools`, `bamana`, and `grammateus`, ONT POD5 example-data policy,
  and synoptikon reporting ambitions.
* Reconciled package metadata and declared currently referenced dependencies so
  package checks can progress beyond initial DESCRIPTION validation.

# floundeR 0.0.5

# floundeR 0.0.4

* inclusion of a basic Rsamtools derived FASTA R6 class (lodestar dependency)
* starting to ensure that methods are also available through documented
  magrittr %>% pipes for legibility of code
* robust Angenieux decorations and demonstration of more comprehensive plots
* integrating Angenieux with ggsave for cleaner Rmarkdown usage

# floundeR 0.0.3

* pkgdown documentation for the github pages
* creation of angenieux prototype; code linking and POC
* creation of the R6 object called FlowCell
* Sequencing summary push ...
* simplified FAST5.is_multi_F5 logic; there were some issues with different
  version of the hdf5 api - tested on Windows and macOS. 

# floundeR 0.0.2

* starting to populate content
* inclusion of FAST5 parsing content, vignette and test data

# floundeR 0.0.1

* tagged a first version

# floundeR 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
