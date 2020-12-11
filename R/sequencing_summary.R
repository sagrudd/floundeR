
#' R6 Class for loading and analysing nanopore sequencing_summary files
#'
#' @description
#' This class aims to simplify the handling and exploration of
#' sequencing_summary files and provides simple methods for accessing
#' information that can be used to assess the performance of a run.
#'
#' @import R6
#' @importFrom utils file_test
#'
#' @export
SequencingSummary <- R6::R6Class(
    "SequencingSummary",
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
        #' sequencing_summary <- system.file("extdata",
        #'     "sequencing_summary.txt.bz2", package="floundeR")
        #' barcodes_summary <- system.file("extdata",
        #'     "barcoding_summary.txt.bz2", package="floundeR")
        #' seqsum <- SequencingSummary$new(
        #'     sequencing_summary_file=sequencing_summary,
        #'     barcoding_summary_file=barcodes_summary)
        initialize = function(
            sequencing_summary_file,
            barcoding_summary_file = NA) {

            # first pass QC - let's ensure that this is a single FILE
            if (is.na(sequencing_summary_file)) {
                stop("Sequencing_summary_file must be defined")
            } else if (!is.character(sequencing_summary_file) ||
                length(sequencing_summary_file) != 1) {
                stop(paste0(
                    "SequencingSummary requires a single [file.path] as input"))
            } else if (!file.exists(sequencing_summary_file)) {
        stop(paste0("path [",sequencing_summary_file,"] does not exist"))
            } else if (!utils::file_test("-f", sequencing_summary_file)) {
                stop(paste0("path [",sequencing_summary_file,
                    "] is a directory - file reqd"))
            }
            self$sequencing_summary_file = sequencing_summary_file

            if (!is.na(barcoding_summary_file)) {
                if (!is.character(barcoding_summary_file) ||
                    length(barcoding_summary_file) != 1) {
                    stop(paste0(
        "barcoding_summary_file requires a single [file.path] as input"))
                } else if (!file.exists(barcoding_summary_file)) {
        stop(paste0("path [",barcoding_summary_file,"] does not exist"))
                } else if (!utils::file_test("-f", barcoding_summary_file)) {
                    stop(paste0("path [",barcoding_summary_file,
                                "] is a directory - file reqd"))
                }
                self$barcoding_summary_file = barcoding_summary_file
            }
            # and import the TSV data (cleanly)
            if (!private$.parse_seqsum()) {
                stop("Failed to import the summary files provided ...")
            }
        }
    ),

    private = list(

        seqsum = NULL,

        .parse_seqsum = function() {
            return(TRUE)
        }

    )
)
