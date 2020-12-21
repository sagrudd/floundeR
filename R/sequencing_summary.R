
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
#'
#' @export
SequencingSummary <- R6::R6Class(
    inherit = FloundeR,
    classname = "SequencingSummary",
    public = list(
        #' @field sequencing_summary_file the file.path
        #' to the query FAST5 file
        sequencing_summary_file = NULL,
        #' @field barcoding_summary_file the file.path
        #' to the query FAST5 file
        barcoding_summary_file = NULL,

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
        #' @return A new `Fast5` object.
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

        sequencingset = function(column="passes_filtering", keys=NULL) {
            if (is.null(private$sequencing_set)) {
                private$sequencing_set <- SequencingSet$new(
                    keycol=column,
                    seqsum=private$seqsum)
            }
            return(private$sequencing_set)
        },

        temporalset = function() {
            if (is.null(private$temporal_set)) {
                private$temporal_set <- TemporalSet$new(
                    seqsum=private$seqsum)
            }
            return(private$temporal_set)
        },

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

        .parse_seqsum = function() {

            mini_table = readr::read_tsv(
                file=self$sequencing_summary_file, n_max=10,
                col_types = readr::cols())
            import_cols <- names(mini_table)[
                na.omit(match(private$select_columns, names(mini_table)))]
            ctypes = private$select_column_types[
                match(import_cols, private$select_columns)]
            colt = paste0(
                "readr::cols_only(",
                paste(paste(
                    import_cols, "=\"", ctypes, "\"",sep=""), collapse=", "),
                ")")

            private$seqsum <- readr::read_tsv(
                file=self$sequencing_summary_file,
                col_types=eval(parse(text=colt)))

            return(TRUE)
        }
    )
)




