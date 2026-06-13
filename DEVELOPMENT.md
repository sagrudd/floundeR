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
