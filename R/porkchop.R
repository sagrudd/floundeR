#' Rank library-preparation kit candidates from read evidence
#'
#' `library_kit_candidates()` calls the curated `../porkchop` Rust library
#' in-process and ranks ONT kit candidates from adapter, primer, barcode, and
#' flank motif evidence. Scores are heuristic evidence scores, not calibrated
#' probabilities. Public default builds keep Porkchop unlinked until that
#' dependency is publicly available; in that mode these wrappers raise typed R
#' conditions.
#'
#' @param reads Character vector of read sequences, or a data frame with
#'   `read_id` and `sequence` columns.
#' @param read_ids Optional character vector of read identifiers. Ignored when
#'   `reads` is a data frame with a `read_id` column.
#'
#' @return A data frame with kit metadata, provenance, support-level,
#'   validation-status, score, `normalized_score`, and `score_kind` columns.
#'
#' @examples
#' if (flounder_rust_capabilities(required = FALSE)$porkchop == "linked") {
#'   library_kit_candidates(c("AATGTACTTCGTTCAGTTACGTATTGCT"))
#' }
#'
#' @export
library_kit_candidates <- function(reads, read_ids = NULL) {
  input <- .flounder_library_reads(reads, read_ids)
  response <- .Call(
    "flounder_library_kit_candidates",
    input$reads,
    input$read_ids,
    PACKAGE = "floundeR"
  )
  .flounder_library_response(response, "Porkchop kit candidate ranking failed.")
}

#' Detect adapter and primer evidence for a selected kit
#'
#' `library_adapter_primer_evidence()` reports adapter and primer motif hits
#' for a selected Porkchop kit. It is read-only QC evidence and does not trim or
#' otherwise transform reads.
#'
#' @param reads Character vector of read sequences, or a data frame with
#'   `read_id` and `sequence` columns.
#' @param kit_id Character scalar. Porkchop kit identifier, such as `LSK114`.
#' @param read_ids Optional character vector of read identifiers. Ignored when
#'   `reads` is a data frame with a `read_id` column.
#'
#' @return A data frame with one row per motif hit.
#'
#' @examples
#' if (flounder_rust_capabilities(required = FALSE)$porkchop == "linked") {
#'   library_adapter_primer_evidence(
#'     reads = c("AATGTACTTCGTTCAGTTACGTATTGCT"),
#'     kit_id = "LSK114"
#'   )
#' }
#'
#' @export
library_adapter_primer_evidence <- function(reads, kit_id, read_ids = NULL) {
  input <- .flounder_library_reads(reads, read_ids)
  kit_id <- .flounder_library_kit_id(kit_id)
  response <- .Call(
    "flounder_library_adapter_primer_evidence",
    input$reads,
    input$read_ids,
    kit_id,
    PACKAGE = "floundeR"
  )
  .flounder_library_motif_response(
    response,
    "Porkchop adapter/primer evidence failed."
  )
}

#' Detect barcode and flank evidence for a selected kit
#'
#' `library_barcode_evidence()` reports barcode and flank motif hits for a
#' selected Porkchop kit. It is read-only QC evidence and does not demultiplex
#' reads or certify barcode assignment.
#'
#' @inheritParams library_adapter_primer_evidence
#'
#' @return A data frame with one row per motif hit.
#'
#' @examples
#' if (flounder_rust_capabilities(required = FALSE)$porkchop == "linked") {
#'   library_barcode_evidence(
#'     reads = c("AAGAAAGTTGTCGGTGTCTTTGTG"),
#'     kit_id = "NBD114.24"
#'   )
#' }
#'
#' @export
library_barcode_evidence <- function(reads, kit_id, read_ids = NULL) {
  input <- .flounder_library_reads(reads, read_ids)
  kit_id <- .flounder_library_kit_id(kit_id)
  response <- .Call(
    "flounder_library_barcode_evidence",
    input$reads,
    input$read_ids,
    kit_id,
    PACKAGE = "floundeR"
  )
  .flounder_library_motif_response(
    response,
    "Porkchop barcode evidence failed."
  )
}

#' Detect cDNA primer-pair evidence for a selected kit
#'
#' `library_cdna_primer_evidence()` reports Porkchop cDNA primer-pair evidence
#' for supported PCS/PCB kits. It preserves Porkchop's class vocabulary and
#' does not imply rescue, demultiplexing, UMI-aware handling, or public-dataset
#' validation.
#'
#' @inheritParams library_adapter_primer_evidence
#'
#' @return A data frame with one row per read and cDNA class/evidence columns.
#'
#' @examples
#' if (flounder_rust_capabilities(required = FALSE)$porkchop == "linked") {
#'   library_cdna_primer_evidence(
#'     reads = c("TTTCTGTTGGTGCTGATATTGCGGCGTCTGCTTGGGTGTTTAACC"),
#'     kit_id = "PCS114"
#'   )
#' }
#'
#' @export
library_cdna_primer_evidence <- function(reads, kit_id, read_ids = NULL) {
  input <- .flounder_library_reads(reads, read_ids)
  kit_id <- .flounder_library_kit_id(kit_id)
  response <- .Call(
    "flounder_library_cdna_primer_evidence",
    input$reads,
    input$read_ids,
    kit_id,
    PACKAGE = "floundeR"
  )
  data <- .flounder_library_response(
    response,
    "Porkchop cDNA primer evidence failed."
  )
  for (column in c(
    "five_prime_name",
    "three_prime_name",
    "workflow_support_note",
    "known_limitations"
  )) {
    data[[column]] <- .flounder_empty_to_na(data[[column]])
  }
  for (column in c(
    "five_prime_start",
    "five_prime_end",
    "three_prime_start",
    "three_prime_end"
  )) {
    data[[column]] <- .flounder_nan_to_na(data[[column]])
  }
  data$primer_hit_count <- as.integer(data$primer_hit_count)
  data
}

.flounder_library_reads <- function(reads, read_ids = NULL) {
  if (is.data.frame(reads)) {
    required <- c("read_id", "sequence")
    missing <- setdiff(required, names(reads))
    if (length(missing) > 0L) {
      .flounder_library_error(
        "`reads` data frame must contain `read_id` and `sequence` columns.",
        category = "argument"
      )
    }
    read_ids <- reads$read_id
    reads <- reads$sequence
  }

  if (!is.character(reads) || anyNA(reads)) {
    .flounder_library_error(
      "`reads` must be a character vector of non-missing sequences.",
      category = "argument"
    )
  }

  if (is.null(read_ids)) {
    read_ids <- sprintf("read_%06d", seq_along(reads))
  }
  if (!is.character(read_ids) || length(read_ids) != length(reads) ||
      anyNA(read_ids)) {
    .flounder_library_error(
      "`read_ids` must be a non-missing character vector matching `reads`.",
      category = "argument"
    )
  }

  list(reads = reads, read_ids = read_ids)
}

.flounder_library_kit_id <- function(kit_id) {
  if (!is.character(kit_id) || length(kit_id) != 1L || is.na(kit_id) ||
      identical(kit_id, "")) {
    .flounder_library_error(
      "`kit_id` must be a non-empty character scalar.",
      category = "argument"
    )
  }
  kit_id
}

.flounder_library_response <- function(response, fallback) {
  if (!is.list(response) || !isTRUE(response$ok)) {
    .flounder_library_error(
      response$error %||% fallback,
      category = response$category %||% "unknown"
    )
  }
  data <- response$data
  for (column in intersect(
    c(
      "provenance_appendix",
      "provenance_notes",
      "source_urls",
      "known_limitations"
    ),
    names(data)
  )) {
    data[[column]] <- .flounder_empty_to_na(data[[column]])
  }
  for (column in intersect(c("introduced_year", "retired_year"), names(data))) {
    data[[column]] <- .flounder_nan_to_na(data[[column]])
  }
  for (column in intersect(c("matched_motifs", "total_hits"), names(data))) {
    data[[column]] <- as.integer(data[[column]])
  }
  data
}

.flounder_library_motif_response <- function(response, fallback) {
  data <- .flounder_library_response(response, fallback)
  data$start <- as.integer(data$start)
  data$end <- as.integer(data$end)
  data$edit_distance <- .flounder_nan_to_na(data$edit_distance)
  data
}

.flounder_library_error <- function(message, category = "unknown", call = NULL) {
  condition <- structure(
    list(message = message, category = category, call = call),
    class = c(
      "floundeR_library_preparation_error",
      "floundeR_error",
      "error",
      "condition"
    )
  )
  stop(condition)
}
