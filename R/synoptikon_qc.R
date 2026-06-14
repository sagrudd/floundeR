#' Build and write Synoptikon QC payloads
#'
#' `as_synoptikon_qc()` assembles existing floundeR QC evidence into the
#' versioned Synoptikon payload contract installed at
#' `inst/schema/synoptikon-qc-payload-v1.schema.json`.
#' `write_synoptikon_qc()` writes that payload as JSON.
#'
#' @param run_summary Optional `qc_run_summary()` one-row tibble.
#' @param report_cards Optional report-card tibble, list of report-card tibbles,
#'   or already wrapped report-card list.
#' @param sequencing_summary,flowcell,barcode,pod5,bam,library_preparation
#'   Optional section evidence. Each may be a data frame/tibble or a named list
#'   of data frames/tibbles.
#' @param input_provenance List of input provenance records.
#' @param run_identity Named list of run, sample, project, flow-cell, kit, and
#'   instrument identifiers.
#' @param governance_context Named list of Synoptikon/Mneion governance
#'   identifiers.
#' @param payload_id Stable payload identifier. When omitted, a timestamped
#'   identifier is generated.
#' @param generated_at_utc Production timestamp as POSIXct or ISO-8601 string.
#' @param limitations List of limitation records.
#' @param intended_use Character vector of Synoptikon intended-use labels.
#' @param compatible_with_grammateus_trusted_reports Logical compatibility flag.
#' @param path Output JSON file path for `write_synoptikon_qc()`.
#' @param pretty Whether to pretty-print JSON.
#' @param ... Passed from `write_synoptikon_qc()` to `as_synoptikon_qc()`.
#'
#' @return `as_synoptikon_qc()` returns a named list compatible with the v1
#'   Synoptikon QC payload schema. `write_synoptikon_qc()` returns `path`
#'   invisibly after writing JSON.
#'
#' @examples
#' summary_file <- flnDr("sequencing_summary.txt.bz2")
#' payload <- as_synoptikon_qc(
#'   run_summary = qc_run_summary(summary_file),
#'   report_cards = qc_report_card(summary_file),
#'   payload_id = "example-qc-payload"
#' )
#'
#' @export
as_synoptikon_qc <- function(
    run_summary = NULL,
    report_cards = NULL,
    sequencing_summary = NULL,
    flowcell = NULL,
    barcode = NULL,
    pod5 = NULL,
    bam = NULL,
    library_preparation = NULL,
    input_provenance = list(),
    run_identity = list(),
    governance_context = list(),
    payload_id = NULL,
    generated_at_utc = Sys.time(),
    limitations = list(),
    intended_use = c("qc_ingestion", "review", "reporting"),
    compatible_with_grammateus_trusted_reports = TRUE) {
  if (is.null(payload_id)) {
    payload_id <- paste0(
      "flounder-qc-",
      format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC", usetz = FALSE))
  }

  run_summary_table <- .synoptikon_table(run_summary)
  run_summary_section <- NULL
  if (!is.null(run_summary_table)) {
    run_summary_section <- list(
      summary = .synoptikon_run_summary_scalars(run_summary_table),
      tables = list(run_summary = run_summary_table))
  }

  sequencing_section <- .synoptikon_section(
    sequencing_summary,
    default = run_summary_section)

  list(
    schema_version = "flounder.synoptikon_qc_payload.v1",
    payload_id = payload_id,
    generated_at_utc = .synoptikon_timestamp(generated_at_utc),
    producer = .synoptikon_producer(),
    run_identity = .synoptikon_named_list(run_identity),
    governance_context = .synoptikon_named_list(governance_context),
    input_provenance = .synoptikon_record_list(input_provenance),
    qc_sections = list(
      sequencing_summary = sequencing_section,
      flowcell = .synoptikon_section(flowcell),
      barcode = .synoptikon_section(barcode),
      pod5 = .synoptikon_section(pod5),
      bam = .synoptikon_section(bam),
      library_preparation = .synoptikon_section(library_preparation)
    ),
    report_cards = .synoptikon_report_cards(report_cards),
    limitations = .synoptikon_record_list(limitations),
    handoff = list(
      consumer = "synoptikon",
      contract_version = "v1",
      intended_use = as.character(intended_use),
      compatible_with_grammateus_trusted_reports =
        isTRUE(compatible_with_grammateus_trusted_reports)
    )
  )
}

#' @rdname as_synoptikon_qc
#' @export
write_synoptikon_qc <- function(path, ..., pretty = TRUE) {
  payload <- as_synoptikon_qc(...)
  json <- jsonlite::toJSON(
    payload,
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    pretty = pretty)
  writeLines(json, path, useBytes = TRUE)
  invisible(path)
}

.synoptikon_producer <- function() {
  list(
    name = "floundeR",
    version = .synoptikon_package_version(),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    rust_capabilities_schema_version = "flounder.rust_capabilities.v1",
    pod5_tools_version = NULL,
    bamana_version = NULL,
    porkchop_version = NULL,
    grammateus_runtime_version = NULL)
}

.synoptikon_package_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("floundeR")),
    error = function(e) {
      description <- file.path(getwd(), "DESCRIPTION")
      if (file.exists(description)) {
        return(read.dcf(description, fields = "Version")[[1]])
      }
      NA_character_
    })
}

.synoptikon_section <- function(x, default = NULL) {
  if (is.null(x) && !is.null(default)) {
    tables <- default$tables
    summary <- default$summary
  } else {
    tables <- .synoptikon_tables(x)
    summary <- .synoptikon_empty_object()
  }

  present <- length(tables) > 0 || length(summary) > 0
  list(
    present = present,
    schema_version = .synoptikon_schema_version(tables),
    status = if (present) "not_checked" else "not_available",
    summary = summary,
    tables = tables,
    artifacts = list())
}

.synoptikon_tables <- function(x) {
  if (is.null(x)) {
    return(.synoptikon_empty_object())
  }

  table <- .synoptikon_table(x)
  if (!is.null(table)) {
    return(list(data = table))
  }

  if (!is.list(x)) {
    stop("Synoptikon section evidence must be a data frame or list.", call. = FALSE)
  }

  out <- list()
  for (name in names(x)) {
    table <- .synoptikon_table(x[[name]])
    if (!is.null(table)) {
      out[[name]] <- table
    }
  }
  if (length(out) == 0) {
    .synoptikon_empty_object()
  } else {
    out
  }
}

.synoptikon_table <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.data.frame(x)) {
    return(NULL)
  }
  lapply(seq_len(nrow(x)), function(i) {
    .synoptikon_named_list(as.list(x[i, , drop = FALSE]))
  })
}

.synoptikon_report_cards <- function(report_cards) {
  if (is.null(report_cards)) {
    return(list())
  }
  if (is.data.frame(report_cards)) {
    return(list(.synoptikon_report_card("report_card", report_cards)))
  }
  if (!is.list(report_cards)) {
    stop("report_cards must be a data frame or list.", call. = FALSE)
  }

  if (!is.null(report_cards$card_id) && !is.null(report_cards$checks)) {
    return(list(report_cards))
  }

  names <- names(report_cards)
  lapply(seq_along(report_cards), function(i) {
    id <- if (!is.null(names) && nzchar(names[[i]])) {
      names[[i]]
    } else {
      paste0("report_card_", i)
    }
    .synoptikon_report_card(id, report_cards[[i]])
  })
}

.synoptikon_report_card <- function(card_id, card) {
  if (!is.data.frame(card)) {
    stop("Each report card must be a data frame.", call. = FALSE)
  }
  rows <- .synoptikon_table(card)
  statuses <- vapply(rows, function(row) row$status %||% NA_character_, character(1))
  list(
    card_id = card_id,
    schema_version = .synoptikon_first(card$schema_version),
    scope = .synoptikon_card_scope(card_id, card),
    status = .synoptikon_worst_status(statuses),
    checks = rows)
}

.synoptikon_card_scope <- function(card_id, card) {
  version <- .synoptikon_first(card$schema_version)
  if (!is.null(version) && grepl("bam", version)) {
    return("bam")
  }
  if (grepl("pod5", card_id)) {
    return("pod5")
  }
  if (grepl("library|porkchop", card_id)) {
    return("library_preparation")
  }
  "run"
}

.synoptikon_worst_status <- function(statuses) {
  statuses <- stats::na.omit(statuses)
  if (length(statuses) == 0) {
    return("not_checked")
  }
  if ("fail" %in% statuses) {
    return("fail")
  }
  if ("warn" %in% statuses) {
    return("warn")
  }
  if (all(statuses == "pass")) {
    return("pass")
  }
  "not_checked"
}

.synoptikon_schema_version <- function(tables) {
  for (table in tables) {
    if (length(table) > 0 && !is.null(table[[1]]$schema_version)) {
      return(table[[1]]$schema_version)
    }
  }
  NULL
}

.synoptikon_run_summary_scalars <- function(rows) {
  if (length(rows) == 0) {
    return(list())
  }
  row <- rows[[1]]
  keep <- intersect(
    names(row),
    c(
      "source_id", "read_count", "passed_read_count", "failed_read_count",
      "pass_fraction", "total_bases", "mean_read_length",
      "n50_read_length", "mean_qscore", "channel_count",
      "run_duration_seconds", "barcode_count", "unclassified_read_count"))
  row[keep]
}

.synoptikon_record_list <- function(x) {
  if (is.null(x)) {
    return(list())
  }
  if (is.data.frame(x)) {
    return(.synoptikon_table(x))
  }
  if (!is.list(x)) {
    stop("Expected a list or data frame.", call. = FALSE)
  }
  if (length(x) == 0) {
    return(list())
  }
  if (is.null(names(x))) {
    return(lapply(x, .synoptikon_named_list))
  }
  list(.synoptikon_named_list(x))
}

.synoptikon_named_list <- function(x) {
  if (is.null(x)) {
    return(.synoptikon_empty_object())
  }
  if (!is.list(x)) {
    stop("Expected a named list.", call. = FALSE)
  }
  if (length(x) == 0) {
    return(.synoptikon_empty_object())
  }
  lapply(x, .synoptikon_scalar)
}

.synoptikon_empty_object <- function() {
  structure(list(), names = character())
}

.synoptikon_scalar <- function(x) {
  if (inherits(x, "POSIXt")) {
    return(.synoptikon_timestamp(x))
  }
  if (inherits(x, "Date")) {
    return(format(x))
  }
  if (length(x) == 0) {
    return(NULL)
  }
  if (length(x) == 1) {
    if (is.factor(x)) {
      return(as.character(x))
    }
    return(x[[1]])
  }
  if (is.factor(x)) {
    return(as.character(x))
  }
  as.vector(x)
}

.synoptikon_first <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NULL)
  }
  x[[1]]
}

.synoptikon_timestamp <- function(x) {
  if (inherits(x, "POSIXt")) {
    return(format(as.POSIXct(x, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"))
  }
  as.character(x)
}

`%||%` <- function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}
