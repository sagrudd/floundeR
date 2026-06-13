# floundeR Development

This file records local development commands that are intentionally kept outside
the built R package.

## Dependency Bootstrap

After FAST5 retirement, `floundeR` no longer depends on `rhdf5`. The current
check surface still needs a mixture of CRAN and Bioconductor packages for the
legacy R code, tests, and vignettes.

Inspect missing dependencies without installing anything:

```sh
Rscript scripts/bootstrap-r-dependencies.R
```

Install missing CRAN and Bioconductor dependencies:

```sh
Rscript scripts/bootstrap-r-dependencies.R --install
```

Fail if any bootstrap dependency is missing:

```sh
Rscript scripts/bootstrap-r-dependencies.R --check
```

The bootstrap script does not install Rust tooling, private Grammateus runtime
assets, ONT POD5 example data, or large external files. Rust source builds need
Cargo and rustc; the development container installs them through rustup.

## Source Install Requirements

Source installs build the embedded Rust extension in `src/rust` through
`src/Makevars` and `src/Makevars.win`. The current minimum Rust toolchain is
Cargo plus rustc with Rust `1.85` or newer, matching `DESCRIPTION` and
`src/rust/Cargo.toml`. Rust-backed APIs are in-process R extension calls, not
CLI wrappers.

The package must remain installable without private Grammateus source or
runtime assets. `../pod5-tools`, `../bamana`, and `../porkchop` are not required
for the current scaffold build; they should become explicit Rust library
dependencies only in the functional binding slices that need them.

### macOS

Install standard build tools and Rust before installing from source:

```sh
xcode-select --install
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
rustc --version
cargo --version
```

For GitHub Actions on macOS, add a Rust toolchain setup step before the R
dependency and package-check steps. The current workflow target is
`macos-latest`; it should ensure `cargo` and `rustc` are on `PATH` before
`R CMD INSTALL` or `R CMD check` runs.

### Linux

Use `rustup` or a system Rust package that provides Rust `1.85` or newer. The
development Docker image uses rustup stable because some distribution packages
lag the Rust edition 2024 floor required by `../pod5-tools`:

```sh
sudo apt-get update
sudo apt-get install -y curl ca-certificates make gcc g++ gfortran pkg-config
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
rustc --version
cargo --version
```

System libraries for the current R dependency surface are declared in
`Dockerfile`. Use that file as the baseline for clean Linux source builds while
the package still depends on native-heavy R packages such as `Rsamtools`,
`ShortRead`, and `magick`.

### Windows

Windows source installs use `src/Makevars.win` and require Rtools plus Cargo
and rustc on `PATH`. The Rust target must match the R/Rtools architecture. This
path is documented but not yet exercised by the project Docker check, so
Windows CI should be added before claiming first-class Windows binary support.

### Docker

The project Docker image is the canonical clean build environment during the
revival:

```sh
docker build -t flounder-dev .
docker run --rm -v "$PWD":/workspace -w /workspace flounder-dev \
  sh scripts/check-r-release-tarball.sh
```

The image installs rustup stable, Cargo, and rustc, then verifies that the
embedded Rust static library links into the R package shared object. It
intentionally does not
include private Grammateus assets, ONT POD5 example data, or large derived
files.

### CI

CI must install or activate Cargo and rustc before R dependency resolution and
package checks. Required CI checks for Rust-bearing changes are:

```sh
cargo fmt --manifest-path=src/rust/Cargo.toml --check
cargo build --manifest-path=src/rust/Cargo.toml --lib --release
Rscript scripts/check-governance-boundaries.R
sh scripts/check-r-release-tarball.sh
```

Credentialed/private CI may add Grammateus runtime checks later, but public CI
must continue to pass without private Grammateus source or runtime downloads.

## Dependency Audit

Write a tab-separated audit of declared package dependencies, their source
channel, installed versions, R/Bioconductor context, and support notes:

```sh
Rscript scripts/audit-r-dependencies.R --output=/tmp/flounder-r-dependencies.tsv
```

The audit exits non-zero when declared dependencies are missing, which makes it
appropriate for CI and container smoke checks.

See `DEPENDENCIES.md` for the current R 4.6/Bioconductor 3.23 container
assessment and package-version snapshot.

## Docker Build

The repository includes a development container that installs the known CRAN
and Bioconductor dependencies into a clean Rocker R image. This keeps dependency
repair reproducible while the package is being revived.

Build the image:

```sh
docker build -t flounder-dev .
```

Run the dependency audit inside the container:

```sh
docker run --rm -v "$PWD":/workspace -w /workspace flounder-dev \
  Rscript scripts/audit-r-dependencies.R --output=/tmp/flounder-r-dependencies.tsv
```

Run release-style package checks inside the container. This builds a source
tarball in a temporary directory and checks that tarball, so diagnostics reflect
the installable package rather than repository-only files such as `.git`,
`.github`, development scripts, and previous check output:

```sh
docker run --rm -v "$PWD":/workspace -w /workspace flounder-dev \
  sh scripts/check-r-release-tarball.sh
```

The container includes rustup stable, Cargo, and rustc so the embedded
`extendr` scaffold can be built during package checks and can later depend on
Rust edition 2024 crates such as `../pod5-tools`. It does not bundle private
Grammateus runtime assets, ONT POD5 example data, or large derived files. Those
remain explicit opt-in layers once the curated Rust interfaces are ready.

## Legacy RMarkdown Vignettes

The historical RMarkdown vignettes are retained under `legacy-vignettes/` as
repo-only reference material while report generation migrates to Grammateus
semantic report contracts. They are intentionally excluded from package builds
and checks because they contain stale plotting side effects, historical
examples, and network-oriented workflows that should not define the rebooted QC
package behavior.

New executable tutorial/report material should be added through the Grammateus
reporting roadmap rather than extending these legacy RMarkdown files.

## Local Checks

Run the repository governance check:

```sh
Rscript scripts/check-governance-boundaries.R
```

Run the package check used during the revival:

```sh
sh scripts/check-r-release-tarball.sh
```

Run focused tests after dependencies are available:

```sh
Rscript -e 'testthat::test_local("tests/testthat")'
```

## Network Tests

Network-dependent tests are disabled by default and must remain safe for
CRAN-like checks. Tests that touch ONT open-data, S3, downloaded POD5 files, or
external services should call `skip_if_no_flounder_network()`.

Run opt-in network tests explicitly:

```sh
FLOUNDER_RUN_NETWORK_TESTS=true Rscript -e 'testthat::test_local("tests/testthat")'
```

Network tests must still use explicit cache directories outside the repository
and must not commit downloaded or derived POD5 files.

## Rust-Backed Tests

Rust-backed tests must skip cleanly when a required compiled feature or runtime
is unavailable. Tests that require compiled POD5, BAM/BGZF/FASTQ, Porkchop, or
Grammateus bindings should call `skip_if_no_flounder_rust()` or a more specific
helper built on `flounder_rust_capabilities()`.

Run Rust-focused tests after installing the package from source:

```sh
Rscript -e 'library(floundeR); testthat::test_file("tests/testthat/test-rust-capabilities.R")'
```
