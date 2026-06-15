#!/usr/bin/env Rscript

truthy <- function(value) {
  identical(tolower(value), "true") || identical(value, "1") ||
    identical(tolower(value), "yes")
}

`%||%` <- function(left, right) {
  if (is.null(left) || length(left) == 0L || is.na(left[[1L]])) {
    return(right)
  }
  left
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

json_value <- function(x, indent = 0L) {
  pad <- paste(rep(" ", indent), collapse = "")
  pad2 <- paste(rep(" ", indent + 2L), collapse = "")
  if (is.data.frame(x)) {
    rows <- lapply(seq_len(nrow(x)), function(i) {
      as.list(x[i, , drop = FALSE])
    })
    return(json_value(rows, indent))
  }
  if (is.list(x) && !is.null(names(x))) {
    keys <- names(x)
    lines <- vapply(seq_along(x), function(i) {
      suffix <- if (i == length(x)) "" else ","
      paste0(
        pad2, "\"", json_escape(keys[[i]]), "\": ",
        json_value(x[[i]], indent + 2L), suffix)
    }, character(1))
    return(paste(c("{", lines, paste0(pad, "}")), collapse = "\n"))
  }
  if (is.list(x)) {
    if (length(x) == 0L) {
      return("[]")
    }
    lines <- vapply(seq_along(x), function(i) {
      suffix <- if (i == length(x)) "" else ","
      paste0(pad2, json_value(x[[i]], indent + 2L), suffix)
    }, character(1))
    return(paste(c("[", lines, paste0(pad, "]")), collapse = "\n"))
  }
  if (length(x) > 1L) {
    return(json_value(as.list(x), indent))
  }
  json_scalar(x)
}

write_json <- function(value, path) {
  writeLines(json_value(value), path)
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

script_source_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) == 0L) {
    return(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  }
  script_path <- sub("^--file=", "", file_arg[[1L]])
  normalizePath(
    file.path(dirname(script_path), ".."),
    winslash = "/",
    mustWork = FALSE)
}

source_root <- script_source_root()

load_flounder_metadata <- function() {
  if (!requireNamespace("floundeR", quietly = TRUE)) {
    source_file <- file.path(source_root, "R", "ont_open_data.R")
    if (!file.exists(source_file)) {
      stop(
        "Install floundeR or run this script from the package source root.",
        call. = FALSE)
    }
    source(source_file, local = .GlobalEnv)
    return(FALSE)
  }
  TRUE
}

flounder_function <- function(name, installed) {
  if (installed) {
    return(getExportedValue("floundeR", name))
  }
  get(name, envir = .GlobalEnv, inherits = FALSE)
}

require_installed_flounder <- function(installed, task) {
  if (!installed) {
    stop(
      task,
      " requires an installed floundeR package with compiled bindings.",
      call. = FALSE)
  }
}

description_field <- function(field) {
  description_path <- file.path(source_root, "DESCRIPTION")
  if (!file.exists(description_path)) {
    return(NA_character_)
  }
  description <- read.dcf(description_path)
  if (!field %in% colnames(description)) {
    return(NA_character_)
  }
  description[1, field]
}

normalise_seqsum_path <- function(path) {
  if (!nzchar(path)) {
    return(NA_character_)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

find_sequence_summary <- function(path) {
  if (!nzchar(path) || !dir.exists(path)) {
    return(NA_character_)
  }
  candidates <- list.files(
    path,
    pattern = "sequencing.*summary|summary.*sequencing",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE)
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0L) {
    return(NA_character_)
  }
  candidates[[which.min(file.info(candidates)$size)]]
}

write_status <- function(path, stage, status, detail) {
  write_tsv(
    data.frame(
      stage = stage,
      status = status,
      detail = detail,
      stringsAsFactors = FALSE),
    path)
}

render_plot_spec <- function(spec, plot_dir, installed) {
  require_installed_flounder(installed, "Rendering real-data QC plots")
  render <- getExportedValue("floundeR", "grammateus_render_plot")
  render(spec, execution = "local_rscript", run_root = plot_dir)
}

build_qc_outputs <- function(seqsum_path, output_root, installed, run_id) {
  require_installed_flounder(installed, "Building real-data QC outputs")
  seqsum <- getExportedValue("floundeR", "SequencingSummary")$new(seqsum_path)
  seqsum_tbl <- seqsum$as_tibble()

  qc_summary <- getExportedValue("floundeR", "qc_run_summary")(
    seqsum_tbl,
    source_id = seqsum_path)
  yield_over_time <- getExportedValue("floundeR", "qc_yield_over_time")(
    seqsum_tbl,
    resolution_minutes = 5)
  read_lengths <- getExportedValue("floundeR", "qc_read_length_distribution")(
    seqsum_tbl,
    bins = 20)
  quality <- getExportedValue("floundeR", "qc_quality_distribution")(
    seqsum_tbl,
    bins = 20)
  channels <- getExportedValue("floundeR", "qc_channel_density")(seqsum_tbl)
  barcodes <- getExportedValue("floundeR", "qc_barcode_composition")(seqsum_tbl)

  table_dir <- file.path(output_root, "qc-tables")
  plot_dir <- file.path(output_root, "plots")
  report_dir <- file.path(output_root, "report")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  tables <- list(
    sequencing_summary = seqsum_tbl,
    qc_summary = qc_summary,
    yield_over_time = yield_over_time,
    read_length_distribution = read_lengths,
    quality_distribution = quality,
    flowcell_density = channels,
    barcode_balance = barcodes
  )
  for (name in names(tables)) {
    write_tsv(tables[[name]], file.path(table_dir, paste0(name, ".tsv")))
  }

  plot_output <- list(width_mm = 140, height_mm = 85, dpi = 150,
                      formats = c("png", "svg"))
  specs <- list(
    yield_over_time = getExportedValue("floundeR", "qc_plot_yield_over_time")(
      yield_over_time, output = plot_output, run_id = run_id),
    quality_distribution =
      getExportedValue("floundeR", "qc_plot_quality_distribution")(
        quality, output = plot_output, run_id = run_id),
    read_length_distribution =
      getExportedValue("floundeR", "qc_plot_read_length_distribution")(
        read_lengths, output = plot_output, run_id = run_id),
    flowcell_density =
      getExportedValue("floundeR", "qc_plot_flowcell_density")(
        channels, output = plot_output, run_id = run_id),
    barcode_balance =
      getExportedValue("floundeR", "qc_plot_barcode_balance")(
        barcodes, output = plot_output, run_id = run_id)
  )
  plot_runs <- lapply(specs, render_plot_spec, plot_dir = plot_dir,
                      installed = installed)
  figure_list <- unlist(lapply(plot_runs, function(run) {
    lapply(run$artifacts, function(artifact) artifact$figure)
  }), recursive = FALSE)

  elements <- getExportedValue("floundeR", "grammateus_qc_report_elements")(
    qc_summary = qc_summary,
    flowcell_density = channels,
    yield_over_time = yield_over_time,
    quality_distribution = quality,
    barcode_balance = barcodes,
    methods = paste(
      "Basecalling is performed outside floundeR. This evidence bundle",
      "ingests a Dorado sequencing-summary file produced by Mnematikon or an",
      "equivalent controlled basecalling workflow, then creates floundeR QC",
      "tables, Grammateus plot artifacts, and a report contract."
    ),
    limitations = paste(
      "This report does not prove basecalling execution by itself. Basecalling",
      "provenance must be supplied by the upstream Mnematikon artifact manifest."
    ),
    provenance = data.frame(
      key = c("sequencing_summary", "basecalling_owner"),
      value = c(normalizePath(seqsum_path, winslash = "/", mustWork = TRUE),
                "Mnematikon or equivalent upstream workflow"),
      stringsAsFactors = FALSE),
    run_id = run_id)

  report <- getExportedValue("floundeR", "qc_report")(
    elements = elements,
    figures = figure_list,
    output_dir = report_dir,
    output = c("html", "pdf"),
    render = "if_available",
    report_id = "report_ont_zymo_pau85136_qc",
    title = "ONT Zymo PAU85136 nanopore QC evidence report",
    run_id = run_id)

  list(
    table_dir = table_dir,
    plot_dir = plot_dir,
    report_dir = report_dir,
    report_contract = report$contract_path,
    report_manifest = report$manifest_path,
    figure_count = length(figure_list),
    table_count = length(tables)
  )
}

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
enabled <- dry_run || truthy(Sys.getenv("FLOUNDER_REAL_DATA_QC", "false"))
if (!enabled) {
  stop(
    "Real-data QC evidence generation is opt-in. Set ",
    "FLOUNDER_REAL_DATA_QC=true or pass --dry-run.",
    call. = FALSE)
}

installed <- load_flounder_metadata()
ont_zymo_pod5_example_objects <- flounder_function(
  "ont_zymo_pod5_example_objects", installed)
objects <- ont_zymo_pod5_example_objects(role = "pass")

output_root <- Sys.getenv(
  "FLOUNDER_REAL_DATA_QC_DIR",
  unset = file.path(
    tools::R_user_dir("floundeR", which = "cache"),
    "real-data-qc-zymo-pau85136"))
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_root)) {
  stop("Could not create output directory: ", output_root, call. = FALSE)
}

metadata_dir <- file.path(output_root, "metadata")
handoff_dir <- file.path(output_root, "mnematikon-handoff")
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(handoff_dir, recursive = TRUE, showWarnings = FALSE)

run_id <- Sys.getenv("FLOUNDER_REAL_DATA_QC_RUN_ID",
                     unset = "ont_zymo_pau85136_pass_8289")
source_path <- Sys.getenv("FLOUNDER_REAL_DATA_QC_POD5", unset = "")
seqsum_path <- Sys.getenv("FLOUNDER_REAL_DATA_QC_SEQUENCE_SUMMARY", unset = "")
basecalled_dir <- Sys.getenv("FLOUNDER_REAL_DATA_QC_BASECALLED_DIR", unset = "")
if (!nzchar(seqsum_path)) {
  seqsum_path <- find_sequence_summary(basecalled_dir)
}

source_metadata_path <- file.path(metadata_dir, "source-metadata.tsv")
write_tsv(objects, source_metadata_path)

source_available <- nzchar(source_path) && file.exists(source_path)
if (source_available && !dry_run) {
  require_installed_flounder(installed, "Recording POD5 raw-data evidence")
  write_tsv(
    getExportedValue("floundeR", "pod5_file_info")(source_path),
    file.path(metadata_dir, "pod5-file-info.tsv"))
  write_tsv(
    getExportedValue("floundeR", "pod5_verify")(source_path),
    file.path(metadata_dir, "pod5-verify.tsv"))
  write_tsv(
    getExportedValue("floundeR", "pod5_manifest")(source_path),
    file.path(metadata_dir, "pod5-manifest.tsv"))
}

handoff <- data.frame(
  field = c(
    "workflow_owner",
    "basecalling_policy",
    "canonical_source",
    "selected_pass_pod5",
    "suggested_run_id",
    "input_pod5_path",
    "basecalled_output_dir",
    "sequencing_summary_path",
    "dgx_host_hint",
    "mnematikon_config_hint"),
  value = c(
    "Mnematikon, not floundeR",
    "floundeR never provides direct basecalling functionality",
    objects$s3_uri[[1]],
    objects$key[[1]],
    run_id,
    normalise_seqsum_path(source_path),
    normalise_seqsum_path(basecalled_dir),
    normalise_seqsum_path(seqsum_path),
    Sys.getenv("FLOUNDER_DGX_HOST", unset = "environment-specific"),
    "../mnematikon"),
  stringsAsFactors = FALSE)
write_tsv(handoff, file.path(handoff_dir, "mnematikon-handoff.tsv"))

analysis_outputs <- NULL
status_path <- file.path(output_root, "analysis-status.tsv")
seqsum_available <- nzchar(seqsum_path) && file.exists(seqsum_path)
if (seqsum_available && !dry_run) {
  analysis_outputs <- build_qc_outputs(
    seqsum_path = seqsum_path,
    output_root = output_root,
    installed = installed,
    run_id = run_id)
  write_status(
    status_path,
    "basecalled_evidence_ingest",
    "complete",
    paste("QC tables, plots, and report contract generated from", seqsum_path))
} else {
  write_status(
    status_path,
    "basecalled_evidence_ingest",
    if (dry_run) "dry-run" else "waiting-for-basecalled-data",
    paste(
      "Provide FLOUNDER_REAL_DATA_QC_SEQUENCE_SUMMARY or",
      "FLOUNDER_REAL_DATA_QC_BASECALLED_DIR after Mnematikon has produced",
      "basecalled outputs. floundeR does not perform basecalling."))
}

manifest <- list(
  schema_version = 1L,
  workflow = "floundeR_real_data_qc_evidence",
  mode = if (dry_run) "dry-run" else "opt-in",
  created_utc = format(as.POSIXct(Sys.time(), tz = "UTC"),
                       "%Y-%m-%dT%H:%M:%SZ"),
  flounder_version = if (installed) {
    as.character(utils::packageVersion("floundeR"))
  } else {
    description_field("Version")
  },
  run_id = run_id,
  source_s3_uri = objects$s3_uri[[1]],
  source_size = objects$size[[1]],
  source_last_modified_utc = objects$last_modified_utc[[1]],
  pod5_path = normalise_seqsum_path(source_path),
  pod5_available = source_available,
  basecalled_dir = normalise_seqsum_path(basecalled_dir),
  sequencing_summary = normalise_seqsum_path(seqsum_path),
  sequencing_summary_available = seqsum_available,
  output_directory = normalizePath(output_root, winslash = "/", mustWork = FALSE),
  basecalling_policy = paste(
    "Basecalling is external to floundeR. Use Mnematikon or another controlled",
    "workflow to generate basecalled data, then pass those artifacts to this",
    "script for QC evidence generation."),
  analysis_outputs = analysis_outputs %||% list()
)
write_json(manifest, file.path(output_root, "workflow-manifest.json"))

readme <- c(
  "# floundeR Real-Data QC Evidence",
  "",
  "This directory is generated outside the package repository.",
  "",
  paste0("- canonical source: `", objects$s3_uri[[1]], "`"),
  "- basecalling owner: Mnematikon or another upstream controlled workflow",
  "- floundeR role: ingest basecalled outputs, produce QC tables, plots, and report contracts",
  "",
  "floundeR does not provide direct basecalling functionality. To build the",
  "full evidence bundle, run Mnematikon against the selected POD5 source and",
  "then rerun this script with either:",
  "",
  "```sh",
  "FLOUNDER_REAL_DATA_QC=true \\",
  "FLOUNDER_REAL_DATA_QC_SEQUENCE_SUMMARY=/path/to/sequencing_summary.tsv \\",
  "Rscript scripts/build-real-data-qc-evidence.R",
  "```",
  "",
  "or set `FLOUNDER_REAL_DATA_QC_BASECALLED_DIR` to a directory containing the",
  "sequencing summary.",
  "",
  "Large POD5, BAM, FASTQ, or derived basecalling artifacts must remain in a",
  "cache or artifact store, not in the floundeR repository."
)
writeLines(readme, file.path(output_root, "README.md"))

message("Real-data QC evidence workflow wrote: ", output_root)
