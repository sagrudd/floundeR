\name{grammateus_plot_spec}
\alias{grammateus_plot_spec}
\alias{grammateus_render_plot}
\title{Build and render Grammateus semantic plot specifications}
\usage{
grammateus_plot_spec(
  plot_id,
  plot_type = c("line", "bar", "stacked_bar", "scatter", "pca"),
  data,
  mappings,
  axes,
  caption,
  theme = "mnemosyne_qc",
  palette = "viridis",
  output = list(width_mm = 160, height_mm = 100, dpi = 300, formats = c("svg",
    "png")),
  bar_value_semantics = NULL,
  source_hash = NULL,
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL
)

grammateus_render_plot(
  plot_spec,
  execution = c("local_rscript", "docker_container"),
  run_root,
  rscript_path = "Rscript",
  docker_path = "docker",
  container_image = NULL,
  container_rscript = "Rscript",
  required_packages = c("ggplot2", "jsonlite", "scales", "svglite", "viridisLite"),
  overwrite_run_dir = TRUE
)
}
\arguments{
\item{plot_id}{Stable lower-snake-case plot identifier. Must start with
\code{plot_}.}

\item{plot_type}{Semantic plot family: \code{line}, \code{bar},
\code{stacked_bar}, \code{scatter}, or \code{pca}.}

\item{data}{Tidy data frame for inline plot data, or a Grammateus data
reference list with \code{kind = "reference"}, \code{path},
\code{source_hash}, and \code{row_count}.}

\item{mappings}{List with required \code{x} and \code{y} fields plus optional
\code{color}, \code{fill}, \code{group}, and \code{label} fields.}

\item{axes}{List with \code{x} and \code{y} axis definitions. Each axis
requires \code{label} and may include \code{unit}, \code{scale}, and
\code{variance_explained}.}

\item{caption}{Scientific plot caption.}

\item{theme}{Named Grammateus theme profile.}

\item{palette}{Named Grammateus palette policy.}

\item{output}{Output definition with \code{width_mm}, \code{height_mm},
\code{dpi}, and \code{formats}.}

\item{bar_value_semantics}{Required for \code{bar} and \code{stacked_bar}
plots. One of \code{counts}, \code{percentages}, \code{rates}, or
\code{normalized_measurements}.}

\item{source_hash}{Optional SHA-256 hash for the source data. Plain
64-character hex digests are accepted and normalised to \code{sha256:<hex>}.}

\item{produced_by}{Producing software or service name.}

\item{producer_version}{Producing software version. When \code{NULL}, the
installed floundeR version is used.}

\item{produced_at_utc}{Production timestamp.}

\item{run_id}{Optional upstream run, workflow, or analysis identifier.}

\item{plot_spec}{A \code{flounder_grammateus_plot_spec} object.}

\item{execution}{Plot backend execution mode. \code{local_rscript} runs an
\code{Rscript} executable directly. \code{docker_container} runs
\code{Rscript} inside a Docker-compatible image with the run directory mounted
at \code{/work}.}

\item{run_root}{Directory where deterministic per-plot run directories are
created.}

\item{rscript_path}{Path to the local \code{Rscript} executable.}

\item{docker_path}{Path to the Docker-compatible CLI executable.}

\item{container_image}{Container image used when \code{execution} is
\code{docker_container}.}

\item{container_rscript}{Rscript command inside the container.}

\item{required_packages}{R packages required by the generated script.}

\item{overwrite_run_dir}{Whether an existing deterministic run directory may
be replaced.}
}
\value{
\code{grammateus_plot_spec()} returns a list with class
\code{flounder_grammateus_plot_spec}. The object follows the Grammateus
\code{ReportPlot} shape and contains plot identity, type, data, mappings,
axes, output settings, optional bar-value semantics, and provenance.

\code{grammateus_render_plot()} returns a list with class
\code{flounder_grammateus_plot_run} containing run paths, hashes, stdout,
stderr, backend session path, and artifact metadata. Artifacts include
\code{flounder_grammateus_figure} wrappers for the generated PNG/SVG files.
}
\description{
\code{grammateus_plot_spec()} creates a backend-neutral Grammateus
\code{ReportPlot}-shaped object from tidy R data. \code{grammateus_render_plot()}
executes the controlled R/ggplot2 backend for that spec, writing deterministic
run inputs, a generated script, backend session metadata, and requested PNG/SVG
artifacts. These helpers keep controlled plot generation available from the
open-source package without requiring private Grammateus runtime assets.
}
