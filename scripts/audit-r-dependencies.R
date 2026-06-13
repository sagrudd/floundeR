#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_arg <- grep("^--output=", args, value = TRUE)
output <- if (length(output_arg) == 0) "" else sub("^--output=", "", output_arg[[1]])

source("scripts/flounder-dependencies.R")
flounder_validate_dependency_plan()

plan <- flounder_dependency_plan()
installed <- as.data.frame(installed.packages(), stringsAsFactors = FALSE)
installed$package <- installed$Package

audit <- merge(plan, installed, by = "package", all.x = TRUE, sort = FALSE)
audit$installed <- !is.na(audit$Version)
audit$installed_version <- ifelse(audit$installed, audit$Version, NA_character_)
audit$built_r_version <- ifelse(audit$installed, audit$Built, NA_character_)
audit$library_path <- ifelse(audit$installed, audit$LibPath, NA_character_)
audit$status <- ifelse(audit$installed, "installed", "missing")

bioc_version <- NA_character_
if (requireNamespace("BiocManager", quietly = TRUE)) {
  bioc_version <- as.character(BiocManager::version())
}

audit$r_version <- paste(R.version$major, R.version$minor, sep = ".")
audit$bioconductor_version <- bioc_version
audit$support_note <- ifelse(
  audit$support_channel == "Bioconductor",
  "Install and update through BiocManager for the active R/Bioconductor release.",
  "Install and update through the configured CRAN mirror."
)

columns <- c(
  "package",
  "declared_in",
  "support_channel",
  "status",
  "installed_version",
  "built_r_version",
  "r_version",
  "bioconductor_version",
  "library_path",
  "support_note"
)
audit <- audit[, columns]
audit <- audit[order(audit$support_channel, audit$package), ]

if (nzchar(output)) {
  write.table(audit, file = output, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Wrote dependency audit:", output, "\n")
} else {
  write.table(audit, file = stdout(), sep = "\t", row.names = FALSE, quote = FALSE)
}

missing <- audit$package[audit$status == "missing"]
if (length(missing) > 0) {
  cat("Missing packages:", paste(missing, collapse = ", "), "\n", file = stderr())
  quit(status = 1)
}
