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
assets, ONT POD5 example data, or large external files.

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

Run package checks inside the container:

```sh
docker run --rm -v "$PWD":/workspace -w /workspace flounder-dev \
  R CMD check --no-manual --no-build-vignettes .
```

The container does not bundle Rust bindings, private Grammateus runtime assets,
ONT POD5 example data, or large derived files. Those remain explicit opt-in
layers once the curated Rust interfaces are ready.

## Local Checks

Run the repository governance check:

```sh
Rscript scripts/check-governance-boundaries.R
```

Run the package check used during the revival:

```sh
R CMD check --no-manual --no-build-vignettes .
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

## Optional Rust Tests

Rust-backed tests must skip cleanly until the in-process Rust extension is
built and available. Tests that require compiled POD5, BAM/BGZF/FASTQ,
Porkchop, or Grammateus bindings should call `skip_if_no_flounder_rust()`.

Run Rust-backed tests explicitly once the compiled bindings exist:

```sh
FLOUNDER_RUST_AVAILABLE=true Rscript -e 'testthat::test_local("tests/testthat")'
```
