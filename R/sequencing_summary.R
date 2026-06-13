
#' R6 Class for loading and analysing nanopore sequencing_summary files
#'
#' @description
#' This class aims to simplify the handling and exploration of
#' sequencing_summary files and provides simple methods for accessing
#' information that can be used to assess the performance of a run.
#'
#' @import R6
#' @importFrom reshape2 acast
#' @importFrom readr read_tsv cols_only
#' @importFrom dplyr count
#' @importFrom magrittr %>%
#'
#' @export
SequencingSummary <- R6::R6Class(
    inherit = FloundeR,
    classname = "SequencingSummary",
    public = list(
        #' @field sequencing_summary_file the file.path
        #' to the query sequencing summary file
        sequencing_summary_file = NULL,
        #' @field barcoding_summary_file the file.path
        #' to the query barcoding summary file
        barcoding_summary_file = NA,

        #' @description
        #' Creates a new SequencingSummary object. This
        #' initialisation method performs other sanity checking
        #' of the defined file(s) to ensure that it is indeed
        #' parseable and creates the required data structures.
        #'
        #' @param sequencing_summary_file The source
        #' sequencing_summary file.
        #' @param barcoding_summary_file The source
        #' barcoding_summary_file file.
        #' @return A new `SequencingSummary` object.
        #'
        #' @examples
        #' sequencing_summary <- flnDr("sequencing_summary.txt.bz2")
        #' barcodes_summary <- flnDr("barcoding_summary.txt.bz2")
        #' seqsum <- SequencingSummary$new(sequencing_summary, barcodes_summary)
        initialize = function(
            sequencing_summary_file,
            barcoding_summary_file = NA) {

            # first pass QC - let's ensure that this is a single FILE
            .check_path(
                "sequencing_summary_file", sequencing_summary_file)
            self$sequencing_summary_file = sequencing_summary_file

            if (!is.na(barcoding_summary_file)) {
                .check_path(
                    "barcoding_summary_file", barcoding_summary_file)
                self$barcoding_summary_file = barcoding_summary_file
            }
            # and import the TSV data (cleanly)
            if (!private$.parse_seqsum()) {
                stop("Failed to import the summary files provided ...")
            }
        },

        #' @description
        #' Export the imported dataset(s) as a tibble
        #'
        #' This object consumes a sequencing summary file (and optionally the
        #' corresponding barcoding_summary file) and creates an object in
        #' memory that can be explored, sliced and filtered. This method dumps
        #' out the in-memory object for further exploration and development.
        #'
        #' @return A tibble representation of the starting dataset
        as_tibble = function() {
            return(private$seqsum)
        }


    ),

    active = list (

        #' @field flowcell
        #' The sequencing summary file contains a collection of data that
        #' includes channel data. The channel data can be overlaid on spatial
        #' representations of flowcell layout to address spatial issues and to
        #' visualise the overall flowcell characteristics. This method creates
        #' a flowcell object that is suitable for these purposes.
        flowcell = function(fc) {

            if (is.null(private$flowcell_object)) {
                private$flowcell_object <- Flowcell$new()
                message("Preparing channel count information")
                channel_counts <- dplyr::count(private$seqsum, channel)
                private$flowcell_object$set_channel_counts(channel_counts)
            }
            return(private$flowcell_object)
        },

        #' @field sequencingset
        #' The `sequencingset` active binding returns a sequencingset object
        #' that is canonically structured around the `passes_filtering` logical
        #' field to allow assessment of sequencing characteristics.
        sequencingset = function(column="passes_filtering", keys=NULL) {
            if (is.null(private$sequencing_set)) {
                private$sequencing_set <- SequencingSet$new(
                    keycol=column,
                    seqsum=private$seqsum)
            }
            return(private$sequencing_set)
        },

        #' @field temporalset
        #' The `temporalset` active binding prepares a temporalset object that
        #' is suitable for the temporal analysis of information within the
        #' sequencing summary file.
        temporalset = function() {
            if (is.null(private$temporal_set)) {
                private$temporal_set <- TemporalSet$new(
                    seqsum=private$seqsum)
            }
            return(private$temporal_set)
        },


        #' @field demultiplex
        #' The `demultiplex` active binding prepared a multiplexset object that
        #' can be used to explore the barcoded content contained within the
        #' sequencing summary file and to access attributes that are related to
        #' these information.
        demultiplex = function() {
            if (is.null(private$multiplex_set)) {
                private$multiplex_set <- MultiplexSet$new(
                    seqsum=private$seqsum,
                    barcoding_summary_file=self$barcoding_summary_file)
            }
            return(private$multiplex_set)
        }

    ),


    private = list(
        flowcell_object = NULL,
        seqsum = NULL,
        sequencing_set = NULL,
        temporal_set = NULL,
        multiplex_set = NULL,

        select_columns = c(
            "read_id", "channel", "start_time", "duration",
            "passes_filtering", "sequence_length_template",
            "mean_qscore_template", "barcode_arrangement"),

        select_column_types = c(
            "c", "i", "d", "d", "l", "d", "d", "c"
        ),

        select_column_aliases = list(
            read_id = c("read_id"),
            channel = c("channel"),
            start_time = c("start_time"),
            duration = c("duration", "template_duration"),
            passes_filtering = c(
                "passes_filtering", "passed_filtering",
                "pass", "passed"),
            sequence_length_template = c(
                "sequence_length_template", "read_length",
                "sequence_length"),
            mean_qscore_template = c(
                "mean_qscore_template", "mean_qscore", "qscore"),
            barcode_arrangement = c("barcode_arrangement", "barcode")
        ),

        optional_columns = c("passes_filtering"),

        .seqsum_col_type = function(type) {
            switch(
                type,
                c = readr::col_character(),
                i = readr::col_integer(),
                d = readr::col_double(),
                l = readr::col_logical(),
                stop(
                    "Unsupported sequencing summary column type: ",
                    type,
                    call. = FALSE))
        },

        .empty_column = function(type, n) {
            switch(
                type,
                c = rep(NA_character_, n),
                i = rep(NA_integer_, n),
                d = rep(NA_real_, n),
                l = rep(NA, n),
                stop(
                    "Unsupported sequencing summary column type: ",
                    type,
                    call. = FALSE))
        },

        .seqsum_col_spec = function(source_cols) {
            ctypes <- private$select_column_types[
                match(source_cols$canonical, private$select_columns)]
            col_specs <- stats::setNames(
                lapply(ctypes, private$.seqsum_col_type),
                source_cols$source)

            do.call(readr::cols_only, col_specs)
        },

        .seqsum_source_columns = function(column_names) {
            rows <- lapply(
                private$select_columns,
                function(canonical) {
                    if (canonical == "read_id" &&
                        is.na(self$barcoding_summary_file)) {
                        return(NULL)
                    }

                    aliases <- private$select_column_aliases[[canonical]]
                    matched <- aliases[aliases %in% column_names]
                    if (length(matched) == 0) {
                        return(NULL)
                    }

                    data.frame(
                        canonical = canonical,
                        source = matched[[1]],
                        stringsAsFactors = FALSE)
                })

            rows <- rows[!vapply(rows, is.null, logical(1))]
            if (length(rows) == 0) {
                return(data.frame(
                    canonical = character(),
                    source = character(),
                    stringsAsFactors = FALSE))
            }

            do.call(rbind, rows)
        },

        .warn_missing_columns = function(missing_cols) {
            if (length(missing_cols) == 0) {
                return(invisible())
            }

            warning(
                "Sequencing summary is missing optional column(s): ",
                paste(missing_cols, collapse = ", "),
                ". Partial QC results will contain NA values for these fields.",
                call. = FALSE)
        },

        .add_missing_columns = function(seqsum, missing_cols) {
            for (column in missing_cols) {
                type <- private$select_column_types[
                    match(column, private$select_columns)]
                seqsum[[column]] <- private$.empty_column(type, nrow(seqsum))
            }

            seqsum
        },

        .parse_seqsum = function() {

            mini_table <- readr::read_tsv(
                file=self$sequencing_summary_file, n_max=10,
                col_types = readr::cols())
            source_cols <- private$.seqsum_source_columns(names(mini_table))
            missing_optional <- setdiff(
                private$optional_columns, source_cols$canonical)
            private$.warn_missing_columns(missing_optional)

            missing_required <- setdiff(
                setdiff(
                    private$select_columns,
                    c(
                        private$optional_columns,
                        "read_id",
                        "barcode_arrangement")),
                source_cols$canonical)
            if (!is.na(self$barcoding_summary_file) &&
                !"read_id" %in% source_cols$canonical) {
                missing_required <- c(missing_required, "read_id")
            }
            if (length(missing_required) > 0) {
                stop(
                    "Sequencing summary is missing required column(s): ",
                    paste(missing_required, collapse = ", "),
                    call. = FALSE)
            }

            private$seqsum <- readr::read_tsv(
                file=self$sequencing_summary_file,
                col_types=private$.seqsum_col_spec(source_cols))
            colnames(private$seqsum) <- source_cols$canonical[
                match(colnames(private$seqsum), source_cols$source)]
            private$seqsum <- private$.add_missing_columns(
                private$seqsum, missing_optional)

            return(TRUE)
        }
    )
)
