#!/usr/bin/env Rscript

truthy <- function(value) {
  identical(tolower(value), "true") || identical(value, "1") ||
    identical(tolower(value), "yes")
}

json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"", "\\\\\"", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  x
}

json_scalar <- function(x) {
  if (length(x) == 0L || is.na(x)) {
    return("null")
  }
  if (is.logical(x)) {
    return(if (isTRUE(x)) "true" else "false")
  }
  if (is.numeric(x)) {
    return(format(x, scientific = FALSE, trim = TRUE))
  }
  paste0("\"", json_escape(as.character(x)), "\"")
}

write_json_object <- function(values, path) {
  keys <- names(values)
  lines <- vapply(seq_along(values), function(i) {
    suffix <- if (i == length(values)) "" else ","
    paste0("  \"", json_escape(keys[[i]]), "\": ", json_scalar(values[[i]]), suffix)
  }, character(1))
  writeLines(c("{", lines, "}"), path)
}

write_tsv <- function(data, path) {
  utils::write.table(
    data,
    file = path,
    sep = "\t",
    row.names = FALSE,
    quote = TRUE,
    na = ""
  )
}

`%||%` <- function(left, right) {
  if (is.null(left) || length(left) == 0L || is.na(left[[1L]])) {
    return(right)
  }
  left
}

description_field <- function(field) {
  if (!file.exists("DESCRIPTION")) {
    return(NA_character_)
  }
  description <- read.dcf("DESCRIPTION")
  if (!field %in% colnames(description)) {
    return(NA_character_)
  }
  description[1, field]
}

pod5_tools_source <- function() {
  cargo_toml <- file.path("src", "rust", "Cargo.toml")
  if (!file.exists(cargo_toml)) {
    return(NA_character_)
  }
  line <- grep("^pod5-tools\\s*=", readLines(cargo_toml, warn = FALSE), value = TRUE)
  if (length(line) == 0L) {
    return(NA_character_)
  }
  line[[1L]]
}

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
enabled <- dry_run || truthy(Sys.getenv("FLOUNDER_DERIVE_POD5_DEMO", "false"))

if (!enabled) {
  stop(
    "Derived POD5 demonstration workflow is opt-in. Set ",
    "FLOUNDER_DERIVE_POD5_DEMO=true or pass --dry-run.",
    call. = FALSE
  )
}

using_installed_package <- requireNamespace("floundeR", quietly = TRUE)
flounder_function <- function(name) {
  if (using_installed_package) {
    return(getExportedValue("floundeR", name))
  }
  get(name, envir = .GlobalEnv, inherits = FALSE)
}

if (!using_installed_package) {
  source_file <- file.path("R", "ont_open_data.R")
  if (!file.exists(source_file)) {
    stop(
      "The floundeR package must be installed, or the script must be run from ",
      "the package source root.",
      call. = FALSE
    )
  }
  source(source_file, local = .GlobalEnv)
}

cache_root <- Sys.getenv(
  "FLOUNDER_DERIVED_POD5_DEMO_DIR",
  unset = file.path(tools::R_user_dir("floundeR", which = "cache"), "derived-pod5-demo")
)
dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(cache_root)) {
  stop("Could not create workflow directory: ", cache_root, call. = FALSE)
}

role <- Sys.getenv("FLOUNDER_DERIVED_POD5_ROLE", unset = "pass")
if (!role %in% c("pass", "fail")) {
  stop("FLOUNDER_DERIVED_POD5_ROLE must be 'pass' or 'fail'.", call. = FALSE)
}

ont_zymo_pod5_example_objects <- flounder_function("ont_zymo_pod5_example_objects")
ont_open_data_fetch <- flounder_function("ont_open_data_fetch")
objects <- ont_zymo_pod5_example_objects(role = role)
source_path <- Sys.getenv("FLOUNDER_DERIVED_POD5_SOURCE", unset = "")
download_source <- truthy(Sys.getenv("FLOUNDER_DERIVE_POD5_DEMO_DOWNLOAD", "false")) &&
  truthy(Sys.getenv("FLOUNDER_RUN_NETWORK_TESTS", "false"))

downloaded <- FALSE
if (!nzchar(source_path) && download_source && !dry_run) {
  fetched <- ont_open_data_fetch(
    key = objects$key[[1]],
    cache_dir = cache_root,
    bucket = objects$bucket[[1]],
    region = objects$region[[1]],
    overwrite = FALSE
  )
  source_path <- fetched$cache_path[[1]]
  downloaded <- isTRUE(fetched$downloaded[[1]])
}

source_available <- nzchar(source_path) && file.exists(source_path)
planned <- FALSE
plan_path <- file.path(cache_root, "subdivide-plan.tsv")
manifest_path <- file.path(cache_root, "pod5-manifest.tsv")

if (source_available && !dry_run) {
  if (!using_installed_package) {
    stop(
      "Local POD5 planning requires an installed floundeR package with Rust support.",
      call. = FALSE
    )
  }
  pod5_subdivide_plan <- flounder_function("pod5_subdivide_plan")
  pod5_manifest <- flounder_function("pod5_manifest")

  plan <- pod5_subdivide_plan(
    source_path,
    strategy = "file-count",
    files_per_chunk = as.integer(Sys.getenv("FLOUNDER_DERIVED_POD5_FILES_PER_CHUNK", "1"))
  )
  write_tsv(plan, plan_path)

  manifest <- pod5_manifest(source_path)
  write_tsv(manifest, manifest_path)
  planned <- TRUE
}

source_metadata_path <- file.path(cache_root, "source-metadata.tsv")
write_tsv(objects, source_metadata_path)

workflow_manifest <- list(
  schema_version = 1L,
  workflow = "floundeR_derived_pod5_demo",
  mode = if (dry_run) "dry-run" else "opt-in",
  created_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"),
  flounder_version = if (using_installed_package) {
    as.character(utils::packageVersion("floundeR"))
  } else {
    description_field("Version")
  },
  pod5_tools_source = pod5_tools_source(),
  source_bucket = objects$bucket[[1]],
  source_key = objects$key[[1]],
  source_s3_uri = objects$s3_uri[[1]],
  source_size = objects$size[[1]],
  source_last_modified_utc = objects$last_modified_utc[[1]],
  source_path = if (nzchar(source_path)) normalizePath(source_path, mustWork = FALSE) else NA_character_,
  source_available = source_available,
  downloaded = downloaded,
  derived_pod5_written = FALSE,
  plan_written = planned,
  manifest_written = planned,
  output_directory = normalizePath(cache_root, mustWork = FALSE),
  policy = paste(
    "This workflow records source metadata and read-only subdivision plans.",
    "It does not write derived POD5 files or commit downloaded data."
  )
)

workflow_manifest_path <- file.path(cache_root, "demo-workflow.json")
write_json_object(workflow_manifest, workflow_manifest_path)

readme_path <- file.path(cache_root, "README.md")
writeLines(c(
  "# floundeR Derived POD5 Demo Workflow",
  "",
  "This directory is generated outside the package repository by an opt-in workflow.",
  "",
  paste0("- source: `", objects$s3_uri[[1]], "`"),
  paste0("- source available locally: `", source_available, "`"),
  paste0("- read-only subdivision plan written: `", planned, "`"),
  "- derived POD5 files written: `FALSE`",
  "",
  "Use `subdivide-plan.tsv` and `pod5-manifest.tsv` as provenance inputs for",
  "documentation, reports, and synoptikon demonstrations. Large source or",
  "derived POD5 files must stay in an explicit cache or artifact store, not in",
  "the floundeR repository."
), readme_path)

message("Derived POD5 demonstration workflow prepared in: ", normalizePath(cache_root, mustWork = FALSE))
message("Workflow manifest: ", workflow_manifest_path)
if (!source_available) {
  message(
    "No local POD5 source was used. Set FLOUNDER_DERIVED_POD5_SOURCE or enable ",
    "FLOUNDER_DERIVE_POD5_DEMO_DOWNLOAD=true with FLOUNDER_RUN_NETWORK_TESTS=true."
  )
}
