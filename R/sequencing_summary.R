
#' R6 Class for loading and analysing nanopore sequencing_summary files
#'
#' @description
#' This class aims to simplify the handling and exploration of
#' sequencing_summary files and provides simple methods for accessing
#' information that can be used to assess the performance of a run.
#'
#' @import R6

#'
#' @export
SequencingSummary <- R6::R6Class("SequencingSummary",
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
                         initialize = function(
                             sequencing_summary_file = NA,
                             barcoding_summary_file = NA) {

                         }
                     )
)
