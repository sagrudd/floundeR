
#' R6 Class for loading and analysing sequence sets
#'
#' @description
#'
#' @importFrom tidyr drop_na
#'
#' @export
SequencingSet <- R6::R6Class(
    inherit = FloundeR,
    classname = "SequencingSet",
    public = list(
        initialize = function(keycol, seqsum=NA) {
            if (is.character(keycol) & tibble::is_tibble(seqsum)) {
                if (! keycol %in% colnames(seqsum)) {
                    stop("keycol must be a seqsum column name")
                }
                private$keycol <- keycol
                private$seqsum <- seqsum
            } else {

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
        },


        read_length_bins = function(normalised=TRUE, cumulative=FALSE, bins=20, outliers=0.025) {
          bins <- self$bin_data(
            private$seqsum$sequence_length_template, bins=bins,
            outliers=outliers)
          private$seqsum["bin"] <- as.numeric(levels(bins))[bins]
          rls <- NULL
          if (normalised) {
            rls <- private$seqsum %>%
              dplyr::group_by(.dots=c(private$keycol, "bin")) %>%
              dplyr::summarize(count=sum(sequence_length_template), .groups="drop")
          } else {
            rls <- private$seqsum %>%
              dplyr::group_by(.dots=c(private$keycol, "bin")) %>%
              dplyr::summarize(count=dplyr::n(), .groups="drop")
          }
          # na may creep into the data ...
          rls <- tidyr::drop_na(rls)

          if (cumulative) {
            rls <- rls %>%
              dplyr::group_by(passes_filtering) %>%
              dplyr::arrange(desc(bin)) %>%
              dplyr::mutate(count=cumsum(count))
          }

          Angenieux$new("2D_count", rls)
        },


        quality_bins = function(bins=20, outliers=0) {
          bins <- self$bin_data(
            private$seqsum$mean_qscore_template, bins=bins,
            outliers=outliers)
          private$seqsum["bin"] <- as.numeric(levels(bins))[bins]
          rls <- private$seqsum %>%
            dplyr::group_by(.dots=c(private$keycol, "bin")) %>%
            dplyr::summarize(count=dplyr::n(), .groups="drop")
          # na may creep into the data ...
          rls <- tidyr::drop_na(rls)

          Angenieux$new("2D_count", rls)
        }
    ),

    active = list (

      enumerate = function(value) {
        if (missing(value)) {
          enumerated_counts <- private$seqsum %>%
            dplyr::group_by(.dots=private$keycol) %>%
            dplyr::summarize(count=dplyr::n())
          Angenieux$new("1D_count", enumerated_counts)
        }
      }

    ),

    private = list(
      keycol = NULL,
      seqsum = NULL
    ),
)
