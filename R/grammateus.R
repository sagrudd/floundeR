#' Build Grammateus figure metadata from image artifacts
#'
#' `grammateus_figure_from_file()` wraps an existing PNG or SVG file as a
#' Grammateus-shaped semantic figure. `grammateus_figure_from_ggplot()` saves a
#' `ggplot2` plot to deterministic PNG and/or SVG artifacts before wrapping the
#' generated files. These helpers do not require the private Grammateus runtime;
#' they prepare the open-source R-side handoff object that later report
#' rendering bindings consume.
#'
#' @param path Path to an existing PNG or SVG file.
#' @param plot A `ggplot2` plot object.
#' @param figure_id Stable lower-snake-case figure identifier. Must start with
#'   `figure_`.
#' @param caption Figure caption.
#' @param alt_text Required accessibility text for HTML output.
#' @param methods_note Optional methods or acquisition note.
#' @param output_dir Directory where plot artifacts are written.
#' @param formats One or more artifact formats: `png`, `svg`.
#' @param width,height Plot dimensions.
#' @param units Unit for `width` and `height`; one of `in`, `cm`, `mm`, or `px`.
#' @param dpi Output resolution used for PNG output and pixel metadata.
#' @param layout_hint Grammateus figure layout hint.
#' @param source_hash Optional SHA-256 hash for the source data that produced
#'   the figure. Plain 64-character hex digests are accepted and normalised to
#'   `sha256:<hex>`.
#' @param source_data Optional source data used to compute `source_hash` when
#'   wrapping a `ggplot2` object. When omitted, the built plot data are used.
#' @param produced_by Producing software or service name.
#' @param producer_version Producing software version.
#' @param produced_at_utc Production timestamp.
#' @param run_id Optional upstream run, workflow, or analysis identifier.
#' @param width_px,height_px Pixel dimensions for an existing file. Required
#'   only when SVG dimensions cannot be inferred.
#' @param transparent_background Whether transparent background preservation is
#'   expected.
#' @param background Background passed to `ggplot2::ggsave()`.
#' @param overwrite Whether existing deterministic artifact paths may be
#'   replaced.
#'
#' @return A list with class `flounder_grammateus_figure`. The object follows
#'   the Grammateus `ReportFigure` shape and contains `figure_id`, `caption`,
#'   `alt_text`, `methods_note`, `layout_hint`, `source`, and `provenance`.
#'
#' @export
grammateus_figure_from_file <- function(
    path,
    figure_id,
    caption,
    alt_text,
    methods_note = NULL,
    layout_hint = c("inline", "full_width", "two_column", "appendix",
                    "panel_group"),
    source_hash = NULL,
    produced_by = "floundeR",
    producer_version = NULL,
    produced_at_utc = Sys.time(),
    run_id = NULL,
    width_px = NULL,
    height_px = NULL,
    transparent_background = FALSE) {
  path <- .grammateus_existing_file(path)
  figure_id <- .grammateus_figure_id(figure_id)
  caption <- .grammateus_required_text(caption, "caption")
  alt_text <- .grammateus_required_text(alt_text, "alt_text")
  layout_hint <- match.arg(layout_hint)
  format <- .grammateus_image_format(path)
  checksum <- .grammateus_file_sha256(path)
  source_hash <- .grammateus_normalize_sha256(
    if (is.null(source_hash)) checksum else source_hash,
    "source_hash"
  )
  if (is.null(producer_version)) {
    producer_version <- .grammateus_default_flounder_version()
  }
  dimensions <- .grammateus_image_dimensions(
    path = path,
    format = format,
    width_px = width_px,
    height_px = height_px
  )

  figure <- list(
    schema_version = "flounder.grammateus_figure.v1",
    figure_id = figure_id,
    caption = caption,
    alt_text = alt_text,
    methods_note = .grammateus_optional_text(methods_note, "methods_note"),
    layout_hint = layout_hint,
    source = list(
      kind = "image",
      format = format,
      path = normalizePath(path, winslash = "/", mustWork = TRUE),
      checksum = checksum,
      width_px = as.integer(dimensions$width_px),
      height_px = as.integer(dimensions$height_px),
      transparent_background = isTRUE(transparent_background)
    ),
    provenance = .grammateus_provenance(
      source_hash = source_hash,
      produced_by = produced_by,
      producer_version = producer_version,
      produced_at_utc = produced_at_utc,
      run_id = run_id
    )
  )
  class(figure) <- c(
    "flounder_grammateus_figure",
    "flounder_grammateus_element",
    "list"
  )
  figure
}

#' @rdname grammateus_figure_from_file
#' @export
grammateus_figure_from_ggplot <- function(
    plot,
    figure_id,
    caption,
    alt_text,
    output_dir,
    formats = c("svg", "png"),
    width = 160,
    height = 100,
    units = c("mm", "in", "cm", "px"),
    dpi = 300,
    methods_note = NULL,
    layout_hint = c("inline", "full_width", "two_column", "appendix",
                    "panel_group"),
    source_hash = NULL,
    source_data = NULL,
    produced_by = "floundeR",
    producer_version = NULL,
    produced_at_utc = Sys.time(),
    run_id = NULL,
    transparent_background = TRUE,
    background = "transparent",
    overwrite = TRUE) {
  if (!inherits(plot, "ggplot")) {
    stop("plot must be a ggplot2 object.", call. = FALSE)
  }
  figure_id <- .grammateus_figure_id(figure_id)
  caption <- .grammateus_required_text(caption, "caption")
  alt_text <- .grammateus_required_text(alt_text, "alt_text")
  methods_note <- .grammateus_optional_text(methods_note, "methods_note")
  layout_hint <- match.arg(layout_hint)
  units <- match.arg(units)
  formats <- .grammateus_plot_formats(formats)
  dimensions <- .grammateus_plot_dimensions(width, height, units, dpi)
  output_dir <- .grammateus_output_dir(output_dir)

  if (is.null(source_hash)) {
    if (is.null(source_data)) {
      source_data <- .grammateus_plot_source_data(plot)
    }
    source_hash <- .grammateus_value_sha256(source_data)
  } else {
    source_hash <- .grammateus_normalize_sha256(source_hash, "source_hash")
  }
  if (is.null(producer_version)) {
    producer_version <- .grammateus_default_flounder_version()
  }

  figures <- lapply(formats, function(format) {
    path <- file.path(output_dir, paste0(figure_id, ".", format))
    if (file.exists(path) && !isTRUE(overwrite)) {
      stop("Artifact already exists: ", path, call. = FALSE)
    }

    if (identical(format, "svg")) {
      if (!requireNamespace("svglite", quietly = TRUE)) {
        stop(
          "The svglite package is required to write SVG Grammateus figures.",
          call. = FALSE
        )
      }
      device <- svglite::svglite
    } else {
      device <- "png"
    }

    ggplot2::ggsave(
      filename = path,
      plot = plot,
      device = device,
      width = dimensions$width_in,
      height = dimensions$height_in,
      units = "in",
      dpi = dpi,
      bg = background
    )

    grammateus_figure_from_file(
      path = path,
      figure_id = figure_id,
      caption = caption,
      alt_text = alt_text,
      methods_note = methods_note,
      layout_hint = layout_hint,
      source_hash = source_hash,
      produced_by = produced_by,
      producer_version = producer_version,
      produced_at_utc = produced_at_utc,
      run_id = run_id,
      width_px = dimensions$width_px,
      height_px = dimensions$height_px,
      transparent_background = transparent_background
    )
  })
  names(figures) <- formats

  if (length(figures) == 1L) {
    return(figures[[1L]])
  }

  bundle <- list(
    schema_version = "flounder.grammateus_figure_bundle.v1",
    figure_id = figure_id,
    caption = caption,
    alt_text = alt_text,
    methods_note = methods_note,
    layout_hint = layout_hint,
    formats = formats,
    figures = figures,
    provenance = figures[[1L]]$provenance
  )
  class(bundle) <- c(
    "flounder_grammateus_figure_bundle",
    "flounder_grammateus_element",
    "list"
  )
  bundle
}

.grammateus_existing_file <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || path == "") {
    stop("path must be a non-empty character scalar.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Figure file does not exist: ", path, call. = FALSE)
  }
  path
}

.grammateus_default_flounder_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("floundeR")),
    error = function(error) "0.0.0-dev"
  )
}

.grammateus_output_dir <- function(output_dir) {
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || output_dir == "") {
    stop("output_dir must be a non-empty character scalar.", call. = FALSE)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Unable to create output_dir: ", output_dir, call. = FALSE)
  }
  normalizePath(output_dir, winslash = "/", mustWork = TRUE)
}

.grammateus_figure_id <- function(figure_id) {
  figure_id <- .grammateus_required_text(figure_id, "figure_id")
  if (!startsWith(figure_id, "figure_")) {
    stop("figure_id must start with `figure_`.", call. = FALSE)
  }
  if (!grepl("^[a-z0-9_]+$", figure_id)) {
    stop("figure_id must use lower snake case.", call. = FALSE)
  }
  figure_id
}

.grammateus_required_text <- function(value, field) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      trimws(value) == "") {
    stop(field, " must be a non-empty character scalar.", call. = FALSE)
  }
  value
}

.grammateus_optional_text <- function(value, field) {
  if (is.null(value)) {
    return(NULL)
  }
  .grammateus_required_text(value, field)
}

.grammateus_image_format <- function(path) {
  format <- tolower(tools::file_ext(path))
  if (!format %in% c("png", "svg")) {
    stop("Figure file must be PNG or SVG.", call. = FALSE)
  }
  format
}

.grammateus_plot_formats <- function(formats) {
  if (!is.character(formats) || length(formats) == 0L || anyNA(formats)) {
    stop("formats must contain at least one format.", call. = FALSE)
  }
  formats <- unique(tolower(formats))
  unsupported <- setdiff(formats, c("png", "svg"))
  if (length(unsupported) > 0L) {
    stop(
      "Unsupported Grammateus figure format(s): ",
      paste(unsupported, collapse = ", "),
      call. = FALSE
    )
  }
  formats
}

.grammateus_plot_dimensions <- function(width, height, units, dpi) {
  for (field in c("width", "height", "dpi")) {
    value <- get(field, inherits = FALSE)
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        value <= 0) {
      stop(field, " must be a positive numeric scalar.", call. = FALSE)
    }
  }
  width_in <- switch(
    units,
    "in" = width,
    "cm" = width / 2.54,
    "mm" = width / 25.4,
    "px" = width / dpi
  )
  height_in <- switch(
    units,
    "in" = height,
    "cm" = height / 2.54,
    "mm" = height / 25.4,
    "px" = height / dpi
  )
  list(
    width_in = width_in,
    height_in = height_in,
    width_px = as.integer(round(width_in * dpi)),
    height_px = as.integer(round(height_in * dpi))
  )
}

.grammateus_file_sha256 <- function(path) {
  digest <- unname(tools::sha256sum(path))
  .grammateus_normalize_sha256(digest, "checksum")
}

.grammateus_value_sha256 <- function(value) {
  path <- tempfile("flounder-grammateus-source-", fileext = ".json")
  on.exit(unlink(path), add = TRUE)
  json <- jsonlite::toJSON(
    value,
    dataframe = "rows",
    null = "null",
    na = "string",
    auto_unbox = TRUE,
    POSIXt = "ISO8601",
    digits = NA
  )
  writeLines(enc2utf8(as.character(json)), path, useBytes = TRUE)
  .grammateus_file_sha256(path)
}

.grammateus_normalize_sha256 <- function(value, field) {
  value <- .grammateus_required_text(value, field)
  value <- tolower(value)
  if (grepl("^[0-9a-f]{64}$", value)) {
    value <- paste0("sha256:", value)
  }
  if (!grepl("^sha256:[0-9a-f]{64}$", value)) {
    stop(
      field,
      " must be a sha256:<hex> digest or 64-character hex digest.",
      call. = FALSE
    )
  }
  value
}

.grammateus_image_dimensions <- function(path, format, width_px = NULL,
                                         height_px = NULL) {
  if (identical(format, "png")) {
    return(.grammateus_png_dimensions(path))
  }

  inferred <- .grammateus_svg_dimensions(path)
  if (!is.null(inferred)) {
    return(inferred)
  }
  if (is.null(width_px) || is.null(height_px)) {
    stop(
      "width_px and height_px are required when SVG dimensions cannot be inferred.",
      call. = FALSE
    )
  }
  .grammateus_positive_pixel_dimensions(width_px, height_px)
}

.grammateus_png_dimensions <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  bytes <- readBin(con, what = "raw", n = 24L)
  signature <- as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))
  if (length(bytes) < 24L || !identical(bytes[1:8], signature)) {
    stop("PNG figure has an invalid header: ", path, call. = FALSE)
  }
  .grammateus_positive_pixel_dimensions(
    .grammateus_uint32_be(bytes[17:20]),
    .grammateus_uint32_be(bytes[21:24])
  )
}

.grammateus_svg_dimensions <- function(path) {
  text <- paste(readLines(path, warn = FALSE, n = 20L), collapse = "\n")
  width <- .grammateus_svg_attr(text, "width")
  height <- .grammateus_svg_attr(text, "height")
  if (is.null(width) || is.null(height)) {
    return(NULL)
  }
  .grammateus_positive_pixel_dimensions(
    .grammateus_svg_length_to_px(width),
    .grammateus_svg_length_to_px(height)
  )
}

.grammateus_svg_attr <- function(text, attr) {
  pattern <- paste0(attr, "\\s*=\\s*['\"]([^'\"]+)['\"]")
  match <- regexec(pattern, text, perl = TRUE)
  value <- regmatches(text, match)[[1]]
  if (length(value) < 2L) {
    return(NULL)
  }
  value[[2L]]
}

.grammateus_svg_length_to_px <- function(value, dpi = 300) {
  numeric <- suppressWarnings(as.numeric(sub(
    "^\\s*([0-9.]+).*$",
    "\\1",
    value,
    perl = TRUE
  )))
  if (is.na(numeric) || numeric <= 0) {
    stop("SVG dimensions must be positive numeric lengths.", call. = FALSE)
  }
  unit <- tolower(sub("^\\s*[0-9.]+\\s*([a-z]*)\\s*$", "\\1", value, perl = TRUE))
  multiplier <- if (identical(unit, "") || identical(unit, "px")) {
    1
  } else {
    switch(
      unit,
      "pt" = dpi / 72,
      "in" = dpi,
      "mm" = dpi / 25.4,
      "cm" = dpi / 2.54,
      1
    )
  }
  as.integer(round(numeric * multiplier))
}

.grammateus_uint32_be <- function(bytes) {
  values <- as.integer(bytes)
  sum(values * c(256^3, 256^2, 256, 1))
}

.grammateus_positive_pixel_dimensions <- function(width_px, height_px) {
  if (!is.numeric(width_px) || length(width_px) != 1L || is.na(width_px) ||
      width_px <= 0 || !is.numeric(height_px) || length(height_px) != 1L ||
      is.na(height_px) || height_px <= 0) {
    stop("Image dimensions must be positive pixel counts.", call. = FALSE)
  }
  list(width_px = as.integer(round(width_px)),
       height_px = as.integer(round(height_px)))
}

.grammateus_plot_source_data <- function(plot) {
  built <- ggplot2::ggplot_build(plot)
  list(
    plot_data = plot$data,
    built_data = built$data
  )
}

.grammateus_provenance <- function(source_hash, produced_by,
                                   producer_version, produced_at_utc,
                                   run_id = NULL) {
  list(
    source_hash = .grammateus_normalize_sha256(source_hash, "source_hash"),
    produced_by = .grammateus_required_text(produced_by, "produced_by"),
    producer_version = .grammateus_required_text(
      producer_version,
      "producer_version"
    ),
    produced_at_utc = .grammateus_utc_timestamp(produced_at_utc),
    run_id = .grammateus_optional_text(run_id, "run_id")
  )
}

.grammateus_utc_timestamp <- function(value) {
  if (inherits(value, "POSIXt")) {
    value <- as.POSIXct(value, tz = "UTC")
    return(format(value, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  }
  .grammateus_required_text(value, "produced_at_utc")
}
