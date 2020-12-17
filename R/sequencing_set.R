
#' R6 Class for loading and analysing sequence sets
#'
#' @description
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
        }
    ),

    active = list (

      set_enumerate = function(value) {
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
