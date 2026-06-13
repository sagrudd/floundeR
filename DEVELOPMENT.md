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
