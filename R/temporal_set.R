
#' R6 Class for analysing sequence sets with accompanying temporal data
#'
#' @description
#'
#' @importFrom tidyr replace_na
#'
#' @export
TemporalSet <- R6::R6Class(
    inherit = FloundeR,
    classname = "TemporalSet",
    public = list(
        initialize = function(seqsum=NA) {
            if (tibble::is_tibble(seqsum)) {
                if (! all(c("start_time", "duration") %in% colnames(seqsum))) {
                    stop("both [start_time] and [duration] columns required")
                }
                private$seqsum <- seqsum
            } else {
                stop("TemporalSet requires tibble as input")
            }
        },

        as_tibble = function() {
            return(private$seqsum)
        }
    ),

    active = list(

    ),

    private = list(
        seqsum = NULL
    )
)
