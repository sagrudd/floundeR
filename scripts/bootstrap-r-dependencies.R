#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
do_install <- "--install" %in% args
do_check <- "--check" %in% args

if (do_install && do_check) {
  stop("Use either --install or --check, not both.", call. = FALSE)
}

source("scripts/flounder-dependencies.R")

cran_packages <- flounder_cran_packages()
bioc_packages <- flounder_bioc_packages()
flounder_validate_dependency_plan()

installed <- rownames(installed.packages())
missing_cran <- setdiff(cran_packages, installed)
missing_bioc <- setdiff(bioc_packages, installed)
missing_all <- c(missing_cran, missing_bioc)

cat("floundeR dependency bootstrap\n")
cat("CRAN packages:", length(cran_packages), "\n")
cat("Bioconductor packages:", length(bioc_packages), "\n")

if (length(missing_all) == 0) {
  cat("All bootstrap dependencies are installed.\n")
  quit(status = 0)
}

cat("Missing CRAN packages:\n")
cat(if (length(missing_cran) == 0) "  none\n" else paste0("  - ", missing_cran, collapse = "\n"), "\n")
cat("Missing Bioconductor packages:\n")
cat(if (length(missing_bioc) == 0) "  none\n" else paste0("  - ", missing_bioc, collapse = "\n"), "\n")

if (do_check) {
  quit(status = 1)
}

if (!do_install) {
  cat("Dry run only. Re-run with --install to install missing packages.\n")
  quit(status = 0)
}

repos <- getOption("repos")
if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}
