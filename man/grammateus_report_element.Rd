\name{grammateus_report_element}
\alias{grammateus_report_element}
\alias{grammateus_qc_report_elements}
\title{Build Grammateus semantic report elements}
\usage{
grammateus_report_element(
  element_id,
  element_type = c("table", "section", "methods", "limitations", "appendix",
    "provenance"),
  title,
  caption = NULL,
  data = NULL,
  body = NULL,
  methods_note = NULL,
  source_hash = NULL,
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL
)

grammateus_qc_report_elements(
  run_metadata = NULL,
  qc_summary = NULL,
  flowcell_density = NULL,
  yield_over_time = NULL,
  quality_distribution = NULL,
  barcode_balance = NULL,
  pod5_integrity = NULL,
  bam_alignment_summary = NULL,
  bam_validation = NULL,
  bam_index = NULL,
  bam_sort = NULL,
  library_preparation = NULL,
  report_card_findings = NULL,
  methods = NULL,
  limitations = NULL,
  appendices = NULL,
  provenance = NULL,
  produced_by = "floundeR",
  producer_version = NULL,
  produced_at_utc = Sys.time(),
  run_id = NULL
)
}
\arguments{
\item{element_id}{Stable lower-snake-case identifier with an element-type prefix such as \code{table_}, \code{section_}, \code{methods_}, \code{limitations_}, \code{appendix_}, or \code{provenance_}.}

\item{element_type}{Semantic element type.}

\item{title}{Human-readable element title.}

\item{caption}{Required caption for table-like report elements.}

\item{data}{Optional data frame or list payload.}

\item{body}{Optional character body for text-like elements.}

\item{methods_note}{Optional methods note.}

\item{source_hash}{Optional SHA-256 hash for the source data. When omitted, a deterministic hash is computed from \code{data} and \code{body}.}

\item{produced_by}{Producing software or service name.}

\item{producer_version}{Producing software version.}

\item{produced_at_utc}{Production timestamp.}

\item{run_id}{Optional upstream run, workflow, or analysis identifier.}

\item{run_metadata, qc_summary, flowcell_density, yield_over_time, quality_distribution, barcode_balance, pod5_integrity, bam_alignment_summary, bam_validation, bam_index, bam_sort, library_preparation, report_card_findings}{Optional QC evidence tables.}

\item{methods}{Optional methods text.}

\item{limitations}{Optional limitations table or text.}

\item{appendices}{Optional named list of appendix tables or text.}

\item{provenance}{Optional provenance table or list.}
}
\value{
A \code{flounder_grammateus_report_element} list, or a \code{flounder_grammateus_report_element_bundle} containing named elements.
}
\description{
\code{grammateus_report_element()} creates a low-level Grammateus-shaped semantic report element from R data. \code{grammateus_qc_report_elements()} builds the standard floundeR QC report element set for run metadata, QC summaries, flowcell/yield/quality/barcode evidence, POD5 integrity, BAM evidence, library-preparation evidence, report-card findings, methods, limitations, appendices, and provenance. These helpers do not render HTML or PDF and do not require private Grammateus runtime assets.
}
