#' Build Grammateus plot specs for nanopore QC report families
#'
#' These helpers convert floundeR QC tables and curated Bamana/POD5 summary
#' objects into backend-neutral Grammateus `ReportPlot` specifications. They do
#' not render plots directly and do not require private Grammateus runtime
#' assets.
#'
#' @param x A floundeR QC table, POD5 evidence table, or `bam_summary()`
#'   result appropriate for the plot family.
#' @param y Metric to plot where the plot family supports alternatives.
#' @param produced_by Producing software or service name.
#' @param producer_version Producing software version.
#' @param produced_at_utc Production timestamp.
#' @param run_id Optional upstream run, workflow, or analysis identifier.
#' @param output Grammateus plot output definition.
#'
#' @return A `flounder_grammateus_plot_spec` object.
#'
#' @export
qc_plot_yield_over_time <- function(
    x,
    y = c("cumulative_bases", "cumulative_read_count"),
    produced_by = "floundeR",
    producer_version = NULL,
    produced_at_utc = Sys.time(),
    run_id = NULL,
    output = NULL) {
  y <- match.arg(y)
  data <- .qc_plot_data_frame(x, "yield_over_time")
  .qc_plot_require(data, c("bin_end_minutes", y, "passes_filtering"))
  data <- .qc_plot_pass_label(data)
  data <- data[, c("bin_end_minutes", y, "pass_state"), drop = FALSE]
  names(data)[names(data) == y] <- "value"
  grammateus_plot_spec(
    plot_id = "plot_yield_over_time",
    plot_type = "line",
    data = data,
    mappings = list(x = "bin_end_minutes", y = "value", color = "pass_state"),
    axes = list(
      x = list(label = "Elapsed time", unit = "minutes"),
      y = list(label = .qc_plot_metric_label(y), unit = .qc_plot_metric_unit(y))
    ),
    caption = "Cumulative sequencing yield over time.",
    output = output %||% .qc_plot_default_output(),
    produced_by = produced_by,
    producer_version = producer_version,
    produced_at_utc = produced_at_utc,
    run_id = run_id
  )
}

#' @rdname qc_plot_yield_over_time
#' @export
qc_plot_quality_distribution <- function(
    x,
    y = c("read_count", "bases"),
    produced_by = "floundeR",
    producer_version = NULL,
    produced_at_utc = Sys.time(),
    run_id = NULL,
    output = NULL) {
  y <- match.arg(y)
  data <- .qc_plot_distribution_data(
    x,
    prefix = "qscore",
    y = y,
    label = "qscore_bin"
  )
  grammateus_plot_spec(
    plot_id = "plot_quality_distribution",
    plot_type = "stacked_bar",
    data = data,
    mappings = list(x = "qscore_bin", y = "value", fill = "pass_state"),
    axes = list(
      x = list(label = "Mean Q-score bin"),
      y = list(label = .qc_plot_metric_label(y), unit = .qc_plot_metric_unit(y))
    ),
    caption = "Distribution of mean read quality.",
    output = output %||% .qc_plot_default_output(),
    bar_value_semantics = "counts",
    produced_by = produced_by,
    producer_version = producer_version,
    produced_at_utc = produced_at_utc,
    run_id = run_id
  )
}

#' @rdname qc_plot_yield_over_time
#' @export
qc_plot_read_length_distribution <- function(
    x,
    y = c("read_count", "bases"),
    produced_by = "floundeR",
    producer_version = NULL,
    produced_at_utc = Sys.time(),
    run_id = NULL,
    output = NULL) {
  y <- match.arg(y)
  data <- .qc_plot_distribution_data(
    x,
    prefix = "read_length",
    y = y,
    label = "read_length_bin"
  )
  grammateus_plot_spec(
    plot_id = "plot_read_length_distribution",
    plot_type = "stacked_bar",
    data = data,
    mappings = list(x = "read_length_bin", y = "value", fill = "pass_state"),
    axes = list(
      x = list(label = "Read length bin", unit = "bases"),
      y = list(label = .qc_plot_metric_label(y), unit = .qc_plot_metric_unit(y))
    ),
    caption = "Distribution of read lengths.",
    output = output %||% .qc_plot_default_output(),
    bar_value_semantics = "counts",
    produced_by = produced_by,
    producer_version = producer_version,
    produced_at_utc = produced_at_utc,
    run_id = run_id
  )
}

#' @rdname qc_plot_yield_over_time
#' @export
qc_plot_flowcell_density <- function(
    x,
    y = c("bases", "read_count", "pass_fraction"),
    produced_by = "floundeR",
    producer_version = NULL,
    produced_at_utc = Sys.time(),
    run_id = NULL,
    output = NULL) {
  y <- match.arg(y)
  data <- .qc_plot_data_frame(x, "flowcell_density")
  .qc_plot_require(data, c("channel", y))
  data <- data[, c("channel", y), drop = FALSE]
  data$channel <- paste0("ch_", data$channel)
  names(data)[names(data) == y] <- "value"
  grammateus_plot_spec(
    plot_id = "plot_flowcell_density",
    plot_type = "bar",
    data = data,
    mappings = list(x = "channel", y = "value"),
    axes = list(
      x = list(label = "Flowcell channel"),
      y = list(label = .qc_plot_metric_label(y), unit = .qc_plot_metric_unit(y))
    ),
    caption = "Flowcell channel density and yield.",
    output = output %||% .qc_plot_default_output(width_mm = 180),
    bar_value_semantics = .qc_plot_bar_semantics(y),
    produced_by = produced_by,
    producer_version = producer_version,
    produced_at_utc = produced_at_utc,
    run_id = run_id
  )
}

#' @rdname qc_plot_yield_over_time
#' @export
qc_plot_barcode_balance <- function(
    x,
    y = c("read_fraction", "bases_fraction", "read_count", "bases"),
    produced_by = "floundeR",
    producer_version = NULL,
    produced_at_utc = Sys.time(),
    run_id = NULL,
    output = NULL) {
  y <- match.arg(y)
  data <- .qc_plot_data_frame(x, "barcode_balance")
  .qc_plot_require(data, c("barcode_arrangement", y))
  data <- data[, c("barcode_arrangement", y), drop = FALSE]
  names(data)[names(data) == y] <- "value"
  grammateus_plot_spec(
    plot_id = "plot_barcode_balance",
    plot_type = "bar",
    data = data,
    mappings = list(x = "barcode_arrangement", y = "value"),
    axes = list(
      x = list(label = "Barcode assignment"),
      y = list(label = .qc_plot_metric_label(y), unit = .qc_plot_metric_unit(y))
    ),
    caption = "Barcode balance by read or base yield.",
    output = output %||% .qc_plot_default_output(width_mm = 160),
    bar_value_semantics = .qc_plot_bar_semantics(y),
    produced_by = produced_by,
    producer_version = producer_version,
    produced_at_utc = produced_at_utc,
    run_id = run_id
  )
}

#' @rdname qc_plot_yield_over_time
#' @export
qc_plot_pod5_integrity <- function(
    x,
    produced_by = "floundeR",
    producer_version = NULL,
    produced_at_utc = Sys.time(),
    run_id = NULL,
    output = NULL) {
  data <- .qc_plot_pod5_integrity_data(x)
  grammateus_plot_spec(
    plot_id = "plot_pod5_integrity",
    plot_type = "bar",
    data = data,
    mappings = list(x = "status", y = "file_count"),
    axes = list(
      x = list(label = "POD5 integrity status"),
      y = list(label = "Files", unit = "count")
    ),
    caption = "POD5 integrity status summary.",
    output = output %||% .qc_plot_default_output(width_mm = 110),
    bar_value_semantics = "counts",
    produced_by = produced_by,
    producer_version = producer_version,
    produced_at_utc = produced_at_utc,
    run_id = run_id
  )
}

#' @rdname qc_plot_yield_over_time
#' @export
qc_plot_bam_mapping_summary <- function(
    x,
    produced_by = "floundeR",
    producer_version = NULL,
    produced_at_utc = Sys.time(),
    run_id = NULL,
    output = NULL) {
  data <- .qc_plot_bam_mapping_data(x)
  grammateus_plot_spec(
    plot_id = "plot_bam_mapping_summary",
    plot_type = "bar",
    data = data,
    mappings = list(x = "category", y = "read_count"),
    axes = list(
      x = list(label = "BAM mapping category"),
      y = list(label = "Reads", unit = "count")
    ),
    caption = "BAM mapping category summary.",
    output = output %||% .qc_plot_default_output(width_mm = 140),
    bar_value_semantics = "counts",
    produced_by = produced_by,
    producer_version = producer_version,
    produced_at_utc = produced_at_utc,
    run_id = run_id
  )
}

#' @rdname qc_plot_yield_over_time
#' @export
qc_plot_bam_mapq_distribution <- function(
    x,
    produced_by = "floundeR",
    producer_version = NULL,
    produced_at_utc = Sys.time(),
    run_id = NULL,
    output = NULL) {
  data <- .qc_plot_bam_mapq_data(x)
  grammateus_plot_spec(
    plot_id = "plot_bam_mapq_distribution",
    plot_type = "bar",
    data = data,
    mappings = list(x = "mapq_label", y = "read_count"),
    axes = list(
      x = list(label = "Mapping quality"),
      y = list(label = "Reads", unit = "count")
    ),
    caption = "BAM mapping-quality distribution.",
    output = output %||% .qc_plot_default_output(width_mm = 150),
    bar_value_semantics = "counts",
    produced_by = produced_by,
    producer_version = producer_version,
    produced_at_utc = produced_at_utc,
    run_id = run_id
  )
}

#' @rdname qc_plot_yield_over_time
#' @export
qc_plot_bam_flag_summary <- function(
    x,
    produced_by = "floundeR",
    producer_version = NULL,
    produced_at_utc = Sys.time(),
    run_id = NULL,
    output = NULL) {
  data <- .qc_plot_bam_flag_data(x)
  grammateus_plot_spec(
    plot_id = "plot_bam_flag_summary",
    plot_type = "bar",
    data = data,
    mappings = list(x = "flag_category", y = "read_count"),
    axes = list(
      x = list(label = "BAM flag category"),
      y = list(label = "Reads", unit = "count")
    ),
    caption = "BAM flag-category summary.",
    output = output %||% .qc_plot_default_output(width_mm = 150),
    bar_value_semantics = "counts",
    produced_by = produced_by,
    producer_version = producer_version,
    produced_at_utc = produced_at_utc,
    run_id = run_id
  )
}

.qc_plot_data_frame <- function(x, label) {
  if (!is.data.frame(x)) {
    stop(label, " plot input must be a data frame.", call. = FALSE)
  }
  as.data.frame(x, stringsAsFactors = FALSE)
}

.qc_plot_require <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    stop(
      "QC plot input is missing column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(data)
}

.qc_plot_default_output <- function(width_mm = 140, height_mm = 90,
                                    dpi = 300) {
  list(width_mm = width_mm, height_mm = height_mm, dpi = dpi,
       formats = c("svg", "png"))
}

.qc_plot_metric_label <- function(metric) {
  labels <- c(
    cumulative_bases = "Cumulative bases",
    cumulative_read_count = "Cumulative reads",
    read_count = "Reads",
    bases = "Bases",
    pass_fraction = "Pass fraction",
    read_fraction = "Read fraction",
    bases_fraction = "Base fraction"
  )
  labels[[metric]] %||% metric
}

.qc_plot_metric_unit <- function(metric) {
  units <- c(
    cumulative_bases = "bases",
    cumulative_read_count = "count",
    read_count = "count",
    bases = "bases",
    pass_fraction = "fraction",
    read_fraction = "fraction",
    bases_fraction = "fraction"
  )
  units[[metric]] %||% "value"
}

.qc_plot_bar_semantics <- function(metric) {
  if (grepl("fraction$", metric)) {
    "rates"
  } else {
    "counts"
  }
}

.qc_plot_pass_label <- function(data) {
  data$pass_state <- ifelse(
    is.na(data$passes_filtering),
    "unknown",
    ifelse(data$passes_filtering %in% TRUE, "passed", "failed")
  )
  data
}

.qc_plot_distribution_data <- function(x, prefix, y, label) {
  data <- .qc_plot_data_frame(x, label)
  start_col <- paste0(prefix, "_bin_start")
  end_col <- paste0(prefix, "_bin_end")
  .qc_plot_require(data, c(start_col, end_col, "passes_filtering", y))
  data <- .qc_plot_pass_label(data)
  out <- data[, c(start_col, end_col, y, "pass_state"), drop = FALSE]
  out[[label]] <- paste0(
    signif(out[[start_col]], 4),
    "-",
    signif(out[[end_col]], 4)
  )
  names(out)[names(out) == y] <- "value"
  out[, c(label, "value", "pass_state"), drop = FALSE]
}

.qc_plot_pod5_integrity_data <- function(x) {
  data <- .qc_plot_data_frame(x, "pod5_integrity")
  if ("verification_status" %in% names(data)) {
    status <- data$verification_status
  } else if ("integrity_status" %in% names(data)) {
    status <- data$integrity_status
  } else if ("overall_status" %in% names(data)) {
    status <- data$overall_status
  } else if ("status" %in% names(data)) {
    status <- data$status
  } else {
    stop(
      "POD5 integrity plot input needs a status or integrity column.",
      call. = FALSE
    )
  }
  status[is.na(status) | status == ""] <- "unknown"
  counts <- as.data.frame(table(status), stringsAsFactors = FALSE)
  names(counts) <- c("status", "file_count")
  counts$file_count <- as.numeric(counts$file_count)
  counts
}

.qc_plot_bam_counts <- function(x) {
  if (inherits(x, "floundeR_bam_summary")) {
    return(.qc_plot_data_frame(x$counts, "bam_counts"))
  }
  .qc_plot_data_frame(x, "bam_counts")
}

.qc_plot_bam_mapping_data <- function(x) {
  counts <- .qc_plot_bam_counts(x)
  fields <- c(
    mapped = "mapped_records",
    unmapped = "unmapped_records",
    primary = "primary_records",
    secondary = "secondary_records",
    supplementary = "supplementary_records"
  )
  .qc_plot_require(counts, unname(fields))
  data.frame(
    category = names(fields),
    read_count = as.numeric(unlist(counts[1, unname(fields), drop = TRUE])),
    stringsAsFactors = FALSE
  )
}

.qc_plot_bam_mapq_data <- function(x) {
  if (inherits(x, "floundeR_bam_summary")) {
    data <- .qc_plot_data_frame(x$mapq_histogram, "bam_mapq")
  } else {
    data <- .qc_plot_data_frame(x, "bam_mapq")
  }
  .qc_plot_require(data, c("mapq", "read_count"))
  data <- data[, c("mapq", "read_count"), drop = FALSE]
  data$mapq_label <- paste0("MAPQ_", data$mapq)
  data$read_count <- as.numeric(data$read_count)
  data[, c("mapq_label", "read_count"), drop = FALSE]
}

.qc_plot_bam_flag_data <- function(x) {
  counts <- .qc_plot_bam_counts(x)
  fields <- c(
    duplicate = "duplicate_records",
    qc_fail = "qc_fail_records",
    paired = "paired_records",
    properly_paired = "properly_paired_records",
    read1 = "read1_records",
    read2 = "read2_records"
  )
  .qc_plot_require(counts, unname(fields))
  data.frame(
    flag_category = names(fields),
    read_count = as.numeric(unlist(counts[1, unname(fields), drop = TRUE])),
    stringsAsFactors = FALSE
  )
}
