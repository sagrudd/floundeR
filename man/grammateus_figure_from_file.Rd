\name{grammateus_figure_from_file}
\alias{grammateus_figure_from_file}
\alias{grammateus_figure_from_ggplot}
\title{Build Grammateus figure metadata from image artifacts}
\usage{
grammateus_figure_from_file(
  path,
  figure_id,
  caption,
  alt_text,
  methods_note = NULL,
  layout_hint = c("inline", "full_width", "two_column", "appendix", "panel_group"),
  source_hash = NULL,
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  width_px = NULL,
  height_px = NULL,
  transparent_background = FALSE
)

grammateus_figure_from_ggplot(
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
  layout_hint = c("inline", "full_width", "two_column", "appendix", "panel_group"),
  source_hash = NULL,
  source_data = NULL,
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL,
  transparent_background = TRUE,
  background = "transparent",
  overwrite = TRUE
)
}
\arguments{
\item{path}{Path to an existing PNG or SVG file.}

\item{figure_id}{Stable lower-snake-case figure identifier. Must start with
\code{figure_}.}

\item{caption}{Figure caption.}

\item{alt_text}{Required accessibility text for HTML output.}

\item{methods_note}{Optional methods or acquisition note.}

\item{layout_hint}{Grammateus figure layout hint.}

\item{source_hash}{Optional SHA-256 hash for the source data that produced
the figure. Plain 64-character hex digests are accepted and normalised to
\code{sha256:<hex>}.}

\item{produced_by}{Producing software or service name.}

\item{producer_version}{Producing software version.}

\item{produced_at_utc}{Production timestamp.}

\item{run_id}{Optional upstream run, workflow, or analysis identifier.}

\item{width_px, height_px}{Pixel dimensions for an existing file. Required
only when SVG dimensions cannot be inferred.}

\item{transparent_background}{Whether transparent background preservation is
expected.}

\item{plot}{A \code{ggplot2} plot object.}

\item{output_dir}{Directory where plot artifacts are written.}

\item{formats}{One or more artifact formats: \code{png}, \code{svg}.}

\item{width, height}{Plot dimensions.}

\item{units}{Unit for \code{width} and \code{height}; one of \code{in},
\code{cm}, \code{mm}, or \code{px}.}

\item{dpi}{Output resolution used for PNG output and pixel metadata.}

\item{source_data}{Optional source data used to compute \code{source_hash}
when wrapping a \code{ggplot2} object. When omitted, the built plot data are
used.}

\item{background}{Background passed to \code{ggplot2::ggsave()}.}

\item{overwrite}{Whether existing deterministic artifact paths may be
replaced.}
}
\value{
A list with class \code{flounder_grammateus_figure}. The object follows the
Grammateus \code{ReportFigure} shape and contains \code{figure_id},
\code{caption}, \code{alt_text}, \code{methods_note}, \code{layout_hint},
\code{source}, and \code{provenance}. When multiple formats are requested from
\code{grammateus_figure_from_ggplot()}, a
\code{flounder_grammateus_figure_bundle} containing one figure per format is
returned.
}
\description{
\code{grammateus_figure_from_file()} wraps an existing PNG or SVG file as a
Grammateus-shaped semantic figure. \code{grammateus_figure_from_ggplot()}
saves a \code{ggplot2} plot to deterministic PNG and/or SVG artifacts before
wrapping the generated files. These helpers do not require the private
Grammateus runtime; they prepare the open-source R-side handoff object that
later report rendering bindings consume.
}
