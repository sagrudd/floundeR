\name{grammateus_runtime_available}
\alias{grammateus_runtime_available}
\alias{grammateus_runtime_version}
\alias{grammateus_runtime_manifest}
\alias{grammateus_runtime_validate}
\alias{grammateus_runtime_install}
\title{Discover and validate a prebuilt Grammateus runtime}
\usage{
grammateus_runtime_available(runtime_root = NULL)

grammateus_runtime_version(runtime_root = NULL)

grammateus_runtime_manifest(runtime_root = NULL)

grammateus_runtime_validate(
  runtime_root = NULL,
  required_capabilities = c("render_report_html", "render_report_pdf")
)

grammateus_runtime_install(
  source_runtime_root,
  cache_root = tools::R_user_dir("floundeR", "cache"),
  overwrite = FALSE
)
}
\arguments{
\item{runtime_root}{Optional explicit runtime root. When omitted, the standard discovery order is used.}

\item{required_capabilities}{Character vector of manifest capabilities that must be \code{TRUE}.}

\item{source_runtime_root}{Existing extracted runtime root to copy into the floundeR-managed cache.}

\item{cache_root}{Optional cache root. Defaults to \code{tools::R_user_dir("floundeR", "cache")}.}

\item{overwrite}{Whether \code{grammateus_runtime_install()} may replace an existing cached runtime directory.}
}
\value{
\code{grammateus_runtime_available()} returns a logical scalar.
\code{grammateus_runtime_version()} returns a character scalar or \code{NA_character_}.
\code{grammateus_runtime_manifest()} returns the parsed manifest list or \code{NULL}.
\code{grammateus_runtime_validate()} returns a validation list with schema version, availability, validity, runtime path, version, platform, capabilities, artifact counts, aggregate artifact hash, and failure records.
\code{grammateus_runtime_install()} returns the validation list for the installed runtime.
}
\description{
These helpers inspect optional private Grammateus runtime bundles without making Grammateus a required dependency of the public floundeR package. Discovery checks \code{GRAMMATEUS_HOME}, then \code{options(floundeR.grammateus_home = "...")}, then the floundeR user cache. No helper downloads private assets during package load or ordinary QC work.
}
