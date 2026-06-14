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

#' Build a Grammateus semantic plot specification
#'
#' `grammateus_plot_spec()` creates a backend-neutral Grammateus `ReportPlot`
#' shaped object from tidy R data. `grammateus_render_plot()` executes the
#' controlled R/ggplot2 backend for that spec, writing deterministic run inputs,
#' a generated script, backend session metadata, and requested PNG/SVG
#' artifacts. These helpers keep controlled plot generation available from the
#' open-source package without requiring private Grammateus runtime assets.
#'
#' @param plot_id Stable lower-snake-case plot identifier. Must start with
#'   `plot_`.
#' @param plot_type Semantic plot family: `line`, `bar`, `stacked_bar`,
#'   `scatter`, or `pca`.
#' @param data Tidy data frame for inline plot data, or a Grammateus data
#'   reference list with `kind = "reference"`, `path`, `source_hash`, and
#'   `row_count`.
#' @param mappings List with required `x` and `y` fields plus optional `color`,
#'   `fill`, `group`, and `label` fields.
#' @param axes List with `x` and `y` axis definitions. Each axis requires
#'   `label` and may include `unit`, `scale`, and `variance_explained`.
#' @param theme Named Grammateus theme profile.
#' @param palette Named Grammateus palette policy.
#' @param output Output definition with `width_mm`, `height_mm`, `dpi`, and
#'   `formats`.
#' @param bar_value_semantics Required for `bar` and `stacked_bar` plots. One
#'   of `counts`, `percentages`, `rates`, or `normalized_measurements`.
#' @param execution Plot backend execution mode. `local_rscript` runs an
#'   `Rscript` executable directly. `docker_container` runs `Rscript` inside a
#'   Docker-compatible image with the run directory mounted at `/work`.
#' @param run_root Directory where deterministic per-plot run directories are
#'   created.
#' @param rscript_path Path to the local `Rscript` executable.
#' @param docker_path Path to the Docker-compatible CLI executable.
#' @param container_image Container image used when `execution` is
#'   `docker_container`.
#' @param container_rscript Rscript command inside the container.
#' @param required_packages R packages required by the generated script.
#' @param overwrite_run_dir Whether an existing deterministic run directory may
#'   be replaced.
#'
#' @return `grammateus_plot_spec()` returns a list with class
#'   `flounder_grammateus_plot_spec`. `grammateus_render_plot()` returns a list
#'   with class `flounder_grammateus_plot_run` containing run paths, hashes,
#'   stdout/stderr, backend session path, and artifact metadata.
#'
#' @export
grammateus_plot_spec <- function(
    plot_id,
    plot_type = c("line", "bar", "stacked_bar", "scatter", "pca"),
    data,
    mappings,
    axes,
    caption,
    theme = "mnemosyne_qc",
    palette = "viridis",
    output = list(width_mm = 160, height_mm = 100, dpi = 300,
                  formats = c("svg", "png")),
    bar_value_semantics = NULL,
    source_hash = NULL,
    produced_by = "floundeR",
    producer_version = NULL,
    produced_at_utc = Sys.time(),
    run_id = NULL) {
  plot_id <- .grammateus_plot_id(plot_id)
  plot_type <- match.arg(plot_type)
  caption <- .grammateus_required_text(caption, "caption")
  mappings <- .grammateus_plot_mappings(mappings)
  axes <- .grammateus_plot_axes(axes, plot_type)
  output <- .grammateus_report_plot_output(output)
  data <- .grammateus_report_plot_data(data, source_hash)
  if (is.null(producer_version)) {
    producer_version <- .grammateus_default_flounder_version()
  }

  if (plot_type %in% c("bar", "stacked_bar") &&
      is.null(bar_value_semantics)) {
    stop(
      "bar and stacked_bar plots must declare bar_value_semantics.",
      call. = FALSE
    )
  }
  if (plot_type %in% c("bar", "stacked_bar")) {
    bar_value_semantics <- .grammateus_bar_value_semantics(
      bar_value_semantics
    )
  } else if (!is.null(bar_value_semantics)) {
    bar_value_semantics <- .grammateus_bar_value_semantics(
      bar_value_semantics
    )
  }
  if (identical(plot_type, "stacked_bar") && is.null(mappings$fill)) {
    stop("stacked_bar plots must declare a fill mapping.", call. = FALSE)
  }

  spec <- list(
    schema_version = "flounder.grammateus_plot.v1",
    plot_id = plot_id,
    plot_type = plot_type,
    caption = caption,
    data = data,
    mappings = mappings,
    axes = axes,
    theme = .grammateus_required_text(theme, "theme"),
    palette = .grammateus_required_text(palette, "palette"),
    output = output,
    provenance = .grammateus_provenance(
      source_hash = .grammateus_plot_data_source_hash(data),
      produced_by = produced_by,
      producer_version = producer_version,
      produced_at_utc = produced_at_utc,
      run_id = run_id
    )
  )
  if (!is.null(bar_value_semantics)) {
    spec$bar_value_semantics <- bar_value_semantics
  }
  .grammateus_validate_plot_spec(spec)
  class(spec) <- c(
    "flounder_grammateus_plot_spec",
    "flounder_grammateus_element",
    "list"
  )
  spec
}

#' @rdname grammateus_plot_spec
#' @export
grammateus_render_plot <- function(
    plot_spec,
    execution = c("local_rscript", "docker_container"),
    run_root,
    rscript_path = "Rscript",
    docker_path = "docker",
    container_image = NULL,
    container_rscript = "Rscript",
    required_packages = c("ggplot2", "jsonlite", "scales", "svglite",
                          "viridisLite"),
    overwrite_run_dir = TRUE) {
  plot_spec <- .grammateus_validate_plot_spec(plot_spec)
  execution <- match.arg(execution)
  run_root <- .grammateus_output_dir(run_root)
  required_packages <- .grammateus_required_package_vector(required_packages)
  run_dir <- file.path(run_root, plot_spec$plot_id)
  if (dir.exists(run_dir)) {
    if (!isTRUE(overwrite_run_dir)) {
      stop("Run directory already exists: ", run_dir, call. = FALSE)
    }
    unlink(run_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(run_dir)) {
    stop("Unable to create run directory: ", run_dir, call. = FALSE)
  }
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)

  spec_path <- file.path(run_dir, "plot_spec.json")
  data_path <- file.path(run_dir, "plot_data.json")
  script_path <- file.path(run_dir, "render_plot.R")
  session_path <- file.path(run_dir, "backend_session.txt")
  .grammateus_write_json(plot_spec, spec_path)
  .grammateus_write_json(plot_spec$data, data_path)
  writeLines(
    .grammateus_r_plot_backend_script(required_packages),
    script_path,
    useBytes = TRUE
  )

  command <- .grammateus_plot_backend_command(
    execution = execution,
    run_dir = run_dir,
    rscript_path = rscript_path,
    docker_path = docker_path,
    container_image = container_image,
    container_rscript = container_rscript
  )
  run_output <- .grammateus_run_backend_command(command, run_dir)
  if (!file.exists(session_path)) {
    .grammateus_abort(
      "R backend did not produce backend_session.txt.",
      "flounder_grammateus_backend_error"
    )
  }

  artifacts <- .grammateus_collect_plot_artifacts(plot_spec, run_dir)
  run <- list(
    schema_version = "flounder.grammateus_plot_run.v1",
    plot_id = plot_spec$plot_id,
    execution = execution,
    run_dir = run_dir,
    script_path = script_path,
    plot_spec_path = spec_path,
    plot_data_path = data_path,
    backend_session_path = session_path,
    plot_spec_sha256 = .grammateus_file_sha256(spec_path),
    plot_data_sha256 = .grammateus_file_sha256(data_path),
    script_sha256 = .grammateus_file_sha256(script_path),
    stdout = run_output$stdout,
    stderr = run_output$stderr,
    artifacts = artifacts,
    provenance = plot_spec$provenance
  )
  class(run) <- c(
    "flounder_grammateus_plot_run",
    "flounder_grammateus_element",
    "list"
  )
  run
}

#' Render Grammateus report elements through the Rust binding
#'
#' `grammateus_render_element()` sends a prepared floundeR Grammateus semantic
#' element to the compiled Rust report-rendering boundary. `format = "html"`
#' returns rendered HTML as a character scalar, while `format = "pdf"` returns
#' PDF bytes as a raw vector. Public floundeR builds do not bundle the private
#' Grammateus runtime; when it is not linked these functions fail with a typed
#' `flounder_grammateus_runtime_unavailable` condition.
#'
#' @param element A Grammateus semantic element prepared by floundeR, currently
#'   a `flounder_grammateus_figure`.
#' @param format Output format: `html` or `pdf`.
#'
#' @return Rendered HTML as a character scalar for `format = "html"`, or PDF
#'   bytes as a raw vector for `format = "pdf"`.
#'
#' @export
grammateus_render_element <- function(element, format = c("html", "pdf")) {
  format <- match.arg(format)
  if (inherits(element, "flounder_grammateus_figure")) {
    if (identical(format, "html")) {
      return(grammateus_render_figure_html(element))
    }
    return(grammateus_render_figure_pdf(element))
  }

  .grammateus_abort(
    paste(
      "Unsupported Grammateus element class:",
      paste(class(element), collapse = ", ")
    ),
    "flounder_grammateus_unsupported_element"
  )
}

#' @rdname grammateus_render_element
#' @export
grammateus_render_figure_html <- function(element) {
  element <- .grammateus_validate_figure(element)
  response <- .Call(
    "flounder_grammateus_render_figure_html",
    element,
    PACKAGE = "floundeR"
  )
  .grammateus_render_response(response, "character")
}

#' @rdname grammateus_render_element
#' @export
grammateus_render_figure_pdf <- function(element) {
  element <- .grammateus_validate_figure(element)
  response <- .Call(
    "flounder_grammateus_render_figure_pdf",
    element,
    PACKAGE = "floundeR"
  )
  .grammateus_render_response(response, "raw")
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

.grammateus_plot_id <- function(plot_id) {
  plot_id <- .grammateus_required_text(plot_id, "plot_id")
  if (!startsWith(plot_id, "plot_")) {
    stop("plot_id must start with `plot_`.", call. = FALSE)
  }
  if (!grepl("^[a-z0-9_]+$", plot_id)) {
    stop("plot_id must use lower snake case.", call. = FALSE)
  }
  plot_id
}

.grammateus_plot_mappings <- function(mappings) {
  if (!is.list(mappings)) {
    stop("mappings must be a list.", call. = FALSE)
  }
  out <- list(
    x = .grammateus_required_text(mappings$x, "x mapping"),
    y = .grammateus_required_text(mappings$y, "y mapping")
  )
  for (field in c("color", "fill", "group", "label")) {
    value <- mappings[[field]]
    if (!is.null(value)) {
      out[[field]] <- .grammateus_required_text(value, paste(field, "mapping"))
    }
  }
  out
}

.grammateus_plot_axes <- function(axes, plot_type) {
  if (!is.list(axes) || !is.list(axes$x) || !is.list(axes$y)) {
    stop("axes must contain x and y axis lists.", call. = FALSE)
  }
  list(
    x = .grammateus_plot_axis(axes$x, "x", plot_type),
    y = .grammateus_plot_axis(axes$y, "y", plot_type)
  )
}

.grammateus_plot_axis <- function(axis, axis_name, plot_type) {
  scale <- axis$scale
  if (is.null(scale)) {
    scale <- "linear"
  }
  scale <- match.arg(scale, c("linear", "log10"))
  out <- list(
    label = .grammateus_required_text(axis$label, paste(axis_name, "axis label")),
    scale = scale
  )
  if (!is.null(axis$unit)) {
    out$unit <- .grammateus_required_text(axis$unit, paste(axis_name, "axis unit"))
  }
  if (!is.null(axis$variance_explained)) {
    if (!is.numeric(axis$variance_explained) ||
        length(axis$variance_explained) != 1L ||
        is.na(axis$variance_explained) ||
        axis$variance_explained < 0 ||
        axis$variance_explained > 1) {
      stop(
        "PCA variance_explained for the ",
        axis_name,
        " axis must be between 0 and 1.",
        call. = FALSE
      )
    }
    out$variance_explained <- axis$variance_explained
  } else if (identical(plot_type, "pca")) {
    stop(
      "PCA plots must declare variance_explained for the ",
      axis_name,
      " axis.",
      call. = FALSE
    )
  }
  out
}

.grammateus_report_plot_output <- function(output) {
  if (!is.list(output)) {
    stop("output must be a list.", call. = FALSE)
  }
  for (field in c("width_mm", "height_mm", "dpi")) {
    value <- output[[field]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        value <= 0) {
      stop("output$", field, " must be a positive numeric scalar.", call. = FALSE)
    }
  }
  if (output$dpi < 72) {
    stop("output$dpi must be at least 72.", call. = FALSE)
  }
  list(
    width_mm = as.integer(round(output$width_mm)),
    height_mm = as.integer(round(output$height_mm)),
    dpi = as.integer(round(output$dpi)),
    formats = .grammateus_plot_formats(output$formats)
  )
}

.grammateus_report_plot_data <- function(data, source_hash = NULL) {
  if (is.data.frame(data)) {
    records <- .grammateus_records_from_data_frame(data)
    if (length(records) == 0L) {
      stop("inline_tidy plot data must contain at least one record.", call. = FALSE)
    }
    out <- list(kind = "inline_tidy", records = records)
    attr(out, "source_hash") <- if (is.null(source_hash)) {
      .grammateus_value_sha256(data)
    } else {
      .grammateus_normalize_sha256(source_hash, "source_hash")
    }
    return(out)
  }

  if (!is.list(data)) {
    stop("data must be a data frame or Grammateus plot data list.", call. = FALSE)
  }
  kind <- .grammateus_required_text(data$kind, "data kind")
  if (identical(kind, "inline_tidy")) {
    if (!is.list(data$records) || length(data$records) == 0L) {
      stop("inline_tidy plot data must contain records.", call. = FALSE)
    }
    attr(data, "source_hash") <- if (is.null(source_hash)) {
      .grammateus_value_sha256(data$records)
    } else {
      .grammateus_normalize_sha256(source_hash, "source_hash")
    }
    return(data)
  }
  if (identical(kind, "reference")) {
    data$path <- .grammateus_required_text(data$path, "data reference path")
    data$source_hash <- .grammateus_normalize_sha256(
      data$source_hash,
      "data reference source_hash"
    )
    if (!is.numeric(data$row_count) || length(data$row_count) != 1L ||
        is.na(data$row_count) || data$row_count <= 0) {
      stop("referenced plot data row_count must be greater than zero.", call. = FALSE)
    }
    data$row_count <- as.integer(round(data$row_count))
    return(data)
  }
  stop("Unsupported plot data kind: ", kind, call. = FALSE)
}

.grammateus_records_from_data_frame <- function(data) {
  rows <- lapply(seq_len(nrow(data)), function(i) {
    row <- lapply(data[i, , drop = FALSE], function(value) {
      value <- value[[1L]]
      if (is.factor(value)) {
        value <- as.character(value)
      }
      if (inherits(value, "POSIXt")) {
        value <- .grammateus_utc_timestamp(value)
      }
      if (is.nan(value)) {
        value <- NA
      }
      if (is.na(value)) {
        return(NULL)
      }
      value
    })
    names(row) <- names(data)
    row
  })
  if (nrow(data) == 0L) {
    list()
  } else {
    rows
  }
}

.grammateus_plot_data_source_hash <- function(data) {
  source_hash <- attr(data, "source_hash", exact = TRUE)
  if (!is.null(source_hash)) {
    return(.grammateus_normalize_sha256(source_hash, "source_hash"))
  }
  if (identical(data$kind, "reference")) {
    return(.grammateus_normalize_sha256(data$source_hash, "source_hash"))
  }
  .grammateus_value_sha256(data)
}

.grammateus_bar_value_semantics <- function(value) {
  choices <- c("counts", "percentages", "rates", "normalized_measurements")
  value <- .grammateus_required_text(value, "bar_value_semantics")
  match.arg(value, choices)
}

.grammateus_validate_plot_spec <- function(spec) {
  if (!is.list(spec)) {
    stop("plot_spec must be a Grammateus plot specification list.", call. = FALSE)
  }
  spec$plot_id <- .grammateus_plot_id(spec$plot_id)
  spec$plot_type <- match.arg(
    .grammateus_required_text(spec$plot_type, "plot_type"),
    c("line", "bar", "stacked_bar", "scatter", "pca")
  )
  spec$caption <- .grammateus_required_text(spec$caption, "caption")
  spec$mappings <- .grammateus_plot_mappings(spec$mappings)
  spec$axes <- .grammateus_plot_axes(spec$axes, spec$plot_type)
  spec$output <- .grammateus_report_plot_output(spec$output)
  spec$theme <- .grammateus_required_text(spec$theme, "theme")
  spec$palette <- .grammateus_required_text(spec$palette, "palette")
  if (spec$plot_type %in% c("bar", "stacked_bar") &&
      is.null(spec$bar_value_semantics)) {
    stop(
      "bar and stacked_bar plots must declare bar_value_semantics.",
      call. = FALSE
    )
  }
  if (!is.null(spec$bar_value_semantics)) {
    spec$bar_value_semantics <- .grammateus_bar_value_semantics(
      spec$bar_value_semantics
    )
  }
  if (identical(spec$plot_type, "stacked_bar") && is.null(spec$mappings$fill)) {
    stop("stacked_bar plots must declare a fill mapping.", call. = FALSE)
  }
  .grammateus_validate_plot_records(spec)
  spec
}

.grammateus_validate_plot_records <- function(spec) {
  if (identical(spec$data$kind, "reference")) {
    invisible(TRUE)
    return()
  }
  if (!identical(spec$data$kind, "inline_tidy") ||
      !is.list(spec$data$records) ||
      length(spec$data$records) == 0L) {
    stop("inline_tidy plot data must contain at least one record.", call. = FALSE)
  }
  required <- unique(unlist(spec$mappings, use.names = FALSE))
  required <- required[!is.na(required) & nzchar(required)]
  for (record_idx in seq_along(spec$data$records)) {
    record <- spec$data$records[[record_idx]]
    for (key in required) {
      if (is.null(record[[key]])) {
        stop(
          "plot record ",
          record_idx - 1L,
          " is missing mapped field `",
          key,
          "`.",
          call. = FALSE
        )
      }
    }
    numeric_keys <- spec$mappings$y
    if (spec$plot_type %in% c("line", "scatter", "pca")) {
      numeric_keys <- c(spec$mappings$x, numeric_keys)
    }
    for (key in unique(numeric_keys)) {
      value <- record[[key]]
      if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
          !is.finite(value)) {
        stop(
          "plot record ",
          record_idx - 1L,
          " mapped field `",
          key,
          "` must be numeric.",
          call. = FALSE
        )
      }
    }
  }
  invisible(TRUE)
}

.grammateus_required_package_vector <- function(required_packages) {
  if (!is.character(required_packages) || length(required_packages) == 0L ||
      anyNA(required_packages) || any(trimws(required_packages) == "")) {
    stop("required_packages must contain at least one package.", call. = FALSE)
  }
  unique(required_packages)
}

.grammateus_write_json <- function(value, path) {
  jsonlite::write_json(
    value,
    path = path,
    dataframe = "rows",
    null = "null",
    na = "null",
    auto_unbox = TRUE,
    POSIXt = "ISO8601",
    pretty = TRUE,
    digits = NA
  )
}

.grammateus_plot_backend_command <- function(execution, run_dir, rscript_path,
                                             docker_path, container_image,
                                             container_rscript) {
  if (identical(execution, "local_rscript")) {
    return(list(command = rscript_path, args = "render_plot.R"))
  }
  container_image <- .grammateus_required_text(
    container_image,
    "container_image"
  )
  list(
    command = docker_path,
    args = c(
      "run", "--rm",
      "-v", paste0(run_dir, ":/work"),
      "-w", "/work",
      container_image,
      container_rscript,
      "render_plot.R"
    )
  )
}

.grammateus_run_backend_command <- function(command, run_dir) {
  stdout_path <- file.path(run_dir, "backend_stdout.txt")
  stderr_path <- file.path(run_dir, "backend_stderr.txt")
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(run_dir)
  status <- suppressWarnings(
    system2(
      command$command,
      args = command$args,
      stdout = stdout_path,
      stderr = stderr_path
    )
  )
  stdout <- if (file.exists(stdout_path)) {
    paste(readLines(stdout_path, warn = FALSE), collapse = "\n")
  } else {
    ""
  }
  stderr <- if (file.exists(stderr_path)) {
    paste(readLines(stderr_path, warn = FALSE), collapse = "\n")
  } else {
    ""
  }
  if (!identical(status, 0L)) {
    .grammateus_abort(
      paste0(
        "R backend failed with status ",
        status,
        if (nzchar(stderr)) paste0(": ", stderr) else "."
      ),
      "flounder_grammateus_backend_error"
    )
  }
  list(stdout = stdout, stderr = stderr)
}

.grammateus_collect_plot_artifacts <- function(plot_spec, run_dir) {
  artifacts <- lapply(plot_spec$output$formats, function(format) {
    path <- file.path(run_dir, paste0("plot.", format))
    if (!file.exists(path)) {
      .grammateus_abort(
        paste0("R backend did not produce artifact: ", path),
        "flounder_grammateus_backend_error"
      )
    }
    file_info <- file.info(path)
    figure <- grammateus_figure_from_file(
      path = path,
      figure_id = sub("^plot_", "figure_", plot_spec$plot_id),
      caption = plot_spec$caption,
      alt_text = paste("Rendered Grammateus plot", plot_spec$plot_id),
      source_hash = plot_spec$provenance$source_hash,
      produced_by = plot_spec$provenance$produced_by,
      producer_version = plot_spec$provenance$producer_version,
      produced_at_utc = plot_spec$provenance$produced_at_utc,
      run_id = plot_spec$provenance$run_id,
      transparent_background = identical(format, "svg")
    )
    list(
      format = format,
      path = normalizePath(path, winslash = "/", mustWork = TRUE),
      sha256 = .grammateus_file_sha256(path),
      byte_len = as.numeric(file_info$size),
      figure = figure
    )
  })
  names(artifacts) <- plot_spec$output$formats
  artifacts
}

.grammateus_validate_figure <- function(figure) {
  if (!inherits(figure, "flounder_grammateus_figure")) {
    stop("element must be a flounder_grammateus_figure.", call. = FALSE)
  }
  figure$figure_id <- .grammateus_figure_id(figure$figure_id)
  figure$caption <- .grammateus_required_text(figure$caption, "caption")
  figure$alt_text <- .grammateus_required_text(figure$alt_text, "alt_text")
  if (!is.list(figure$source)) {
    stop("figure source must be a list.", call. = FALSE)
  }
  figure$source$kind <- .grammateus_required_text(
    figure$source$kind,
    "figure source kind"
  )
  if (!identical(figure$source$kind, "image")) {
    stop("figure source kind must be `image`.", call. = FALSE)
  }
  figure$source$format <- match.arg(figure$source$format, c("png", "svg"))
  figure$source$path <- .grammateus_existing_file(figure$source$path)
  figure$source$checksum <- .grammateus_normalize_sha256(
    figure$source$checksum,
    "figure source checksum"
  )
  if (!is.list(figure$provenance)) {
    stop("figure provenance must be a list.", call. = FALSE)
  }
  figure$provenance$source_hash <- .grammateus_normalize_sha256(
    figure$provenance$source_hash,
    "figure provenance source_hash"
  )
  figure
}

.grammateus_render_response <- function(response, expected_type) {
  if (!is.list(response) || is.null(response$ok)) {
    .grammateus_abort(
      "Grammateus Rust renderer returned an invalid response envelope.",
      "flounder_grammateus_render_error"
    )
  }
  if (isTRUE(response$ok)) {
    data <- response$data
    if (identical(expected_type, "character")) {
      if (!is.character(data) || length(data) != 1L || is.na(data)) {
        .grammateus_abort(
          "Grammateus Rust renderer did not return HTML text.",
          "flounder_grammateus_render_error"
        )
      }
      return(data)
    }
    if (!is.raw(data)) {
      .grammateus_abort(
        "Grammateus Rust renderer did not return PDF bytes.",
        "flounder_grammateus_render_error"
      )
    }
    return(data)
  }

  category <- response$category
  if (!is.character(category) || length(category) != 1L || is.na(category)) {
    category <- "render"
  }
  class <- if (identical(category, "runtime_unavailable")) {
    "flounder_grammateus_runtime_unavailable"
  } else {
    "flounder_grammateus_render_error"
  }
  message <- response$error
  if (!is.character(message) || length(message) != 1L || is.na(message)) {
    message <- "Grammateus Rust renderer failed."
  }
  .grammateus_abort(message, class)
}

.grammateus_r_plot_backend_script <- function(required_packages) {
  packages <- paste(
    vapply(required_packages, function(package) {
      paste0('"', gsub('"', '\\"', package, fixed = TRUE), '"')
    }, character(1)),
    collapse = ", "
  )
  paste0(
    "# floundeR Grammateus controlled R/ggplot2 backend template v1\n",
    "required_packages <- c(", packages, ")\n",
    "missing_packages <- required_packages[!vapply(required_packages, ",
    "requireNamespace, logical(1), quietly = TRUE)]\n",
    "if (length(missing_packages) > 0) {\n",
    "  stop(paste('missing required R packages:', ",
    "paste(missing_packages, collapse = ', ')), call. = FALSE)\n",
    "}\n\n",
    "plot_spec <- jsonlite::fromJSON('plot_spec.json', simplifyVector = FALSE)\n",
    "plot_data <- jsonlite::fromJSON('plot_data.json', simplifyDataFrame = TRUE)\n",
    "width_in <- plot_spec$output$width_mm / 25.4\n",
    "height_in <- plot_spec$output$height_mm / 25.4\n",
    "dpi <- plot_spec$output$dpi\n\n",
    "if (identical(plot_data$kind, 'inline_tidy')) {\n",
    "  data <- as.data.frame(plot_data$records, stringsAsFactors = FALSE)\n",
    "} else if (identical(plot_data$kind, 'reference')) {\n",
    "  data <- utils::read.csv(plot_data$path, stringsAsFactors = FALSE)\n",
    "} else {\n",
    "  stop(paste('unsupported plot data kind:', plot_data$kind), call. = FALSE)\n",
    "}\n\n",
    "mapping <- plot_spec$mappings\n",
    "aes_args <- list(x = mapping$x, y = mapping$y)\n",
    "if (!is.null(mapping$color)) aes_args$colour <- mapping$color\n",
    "if (!is.null(mapping$fill)) aes_args$fill <- mapping$fill\n",
    "if (!is.null(mapping$group)) aes_args$group <- mapping$group\n",
    "if (!is.null(mapping$label)) aes_args$label <- mapping$label\n\n",
    "axis_label <- function(axis) {\n",
    "  label <- axis$label\n",
    "  if (!is.null(axis$variance_explained)) {\n",
    "    label <- paste0(label, ' (', sprintf('%.1f', ",
    "axis$variance_explained * 100), '%)')\n",
    "  } else if (!is.null(axis$unit)) {\n",
    "    label <- paste0(label, ' (', axis$unit, ')')\n",
    "  }\n",
    "  label\n",
    "}\n\n",
    "plot <- ggplot2::ggplot(data, do.call(ggplot2::aes_string, aes_args))\n",
    "if (identical(plot_spec$plot_type, 'line')) {\n",
    "  plot <- plot + ggplot2::geom_line(linewidth = 0.6) + ",
    "ggplot2::geom_point(size = 1.8)\n",
    "} else if (identical(plot_spec$plot_type, 'bar')) {\n",
    "  plot <- plot + ggplot2::geom_col(width = 0.72)\n",
    "} else if (identical(plot_spec$plot_type, 'stacked_bar')) {\n",
    "  plot <- plot + ggplot2::geom_col(width = 0.72, position = 'stack')\n",
    "} else if (identical(plot_spec$plot_type, 'scatter')) {\n",
    "  plot <- plot + ggplot2::geom_point(size = 2.4, alpha = 0.85)\n",
    "} else if (identical(plot_spec$plot_type, 'pca')) {\n",
    "  plot <- plot + ggplot2::geom_point(size = 2.6, alpha = 0.9)\n",
    "  if (!is.null(mapping$label)) {\n",
    "    plot <- plot + ggplot2::geom_text(vjust = -0.8, size = 2.6, ",
    "show.legend = FALSE)\n",
    "  }\n",
    "} else {\n",
    "  stop(paste('unsupported plot type:', plot_spec$plot_type), call. = FALSE)\n",
    "}\n\n",
    "if (!is.null(mapping$colour) || !is.null(mapping$color)) {\n",
    "  color_key <- if (!is.null(mapping$colour)) mapping$colour else mapping$color\n",
    "  color_levels <- length(unique(data[[color_key]]))\n",
    "  plot <- plot + ggplot2::scale_color_manual(values = ",
    "viridisLite::viridis(color_levels))\n",
    "}\n",
    "if (!is.null(mapping$fill)) {\n",
    "  fill_levels <- length(unique(data[[mapping$fill]]))\n",
    "  plot <- plot + ggplot2::scale_fill_manual(values = ",
    "viridisLite::viridis(fill_levels))\n",
    "}\n\n",
    "plot <- plot +\n",
    "  ggplot2::labs(\n",
    "    title = plot_spec$caption,\n",
    "    x = axis_label(plot_spec$axes$x),\n",
    "    y = axis_label(plot_spec$axes$y),\n",
    "    color = mapping$color,\n",
    "    fill = mapping$fill\n",
    "  ) +\n",
    "  ggplot2::theme_minimal(base_size = 10) +\n",
    "  ggplot2::theme(\n",
    "    plot.title = ggplot2::element_text(face = 'bold', size = 11),\n",
    "    panel.grid.minor = ggplot2::element_blank(),\n",
    "    legend.position = 'right'\n",
    "  )\n\n",
    "if ('png' %in% plot_spec$output$formats) {\n",
    "  grDevices::png('plot.png', width = width_in, height = height_in, ",
    "units = 'in', res = dpi)\n",
    "  print(plot)\n",
    "  grDevices::dev.off()\n",
    "}\n\n",
    "if ('svg' %in% plot_spec$output$formats) {\n",
    "  svglite::svglite('plot.svg', width = width_in, height = height_in)\n",
    "  print(plot)\n",
    "  grDevices::dev.off()\n",
    "}\n\n",
    "writeLines(capture.output(sessionInfo()), 'backend_session.txt')\n",
    "invisible(plot_data)\n"
  )
}

.grammateus_abort <- function(message, class) {
  stop(
    structure(
      list(message = message, call = NULL),
      class = c(class, "flounder_grammateus_error", "error", "condition")
    )
  )
}
