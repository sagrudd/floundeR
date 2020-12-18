
#' R6 Class for analysing sequence sets with accompanying temporal data
#'
#' @description
#'
#' @importFrom tidyr drop_na
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
        },

        runtime = function() {
            runtime <- max(private$seqsum$start_time) / 60 / 60
            expectedRuntimes <- c(4,8,12,24,36,48,64,72,96)
            temporaldistance <- sqrt((expectedRuntimes - runtime)^2)
            rruntime <- expectedRuntimes[[which(
                temporaldistance==min(temporaldistance))]]
            return(rruntime)
        },

        read_yield = function(resolution=15, cumulative=TRUE) {
            bins <- self$bin_data(
                private$seqsum$start_time/60/60,
                bins=self$runtime()*60/resolution,
                qmax=self$runtime())
            # bins is thus in hours with resolution in minutes ...
            private$seqsum["bin"] <- as.numeric(levels(bins))[bins]
            rls <- private$seqsum %>%
                dplyr::group_by(.dots=c("passes_filtering", "bin")) %>%
                dplyr::summarize(count=dplyr::n(), .groups="drop")
            rls <- tidyr::drop_na(rls)

            if (cumulative) {
                rls <- rls %>% group_by(passes_filtering) %>% arrange(bin) %>% mutate(count=cumsum(count))
            }

            Angenieux$new("2D_count", rls)
        },

        base_yield = function(resolution=15, cumulative=TRUE) {
            bins <- self$bin_data(
                private$seqsum$start_time/60/60,
                bins=self$runtime()*60/resolution,
                qmax=self$runtime())
            # bins is thus in hours with resolution in minutes ...
            private$seqsum["bin"] <- as.numeric(levels(bins))[bins]
            rls <- private$seqsum %>%
                dplyr::group_by(.dots=c("passes_filtering", "bin")) %>%
                dplyr::summarize(count=sum(sequence_length_template), .groups="drop")
            rls <- tidyr::drop_na(rls)

            if (cumulative) {
                rls <- rls %>% group_by(passes_filtering) %>% arrange(bin) %>% mutate(count=cumsum(count))
            }

            Angenieux$new("2D_count", rls)
        }


    ),

    active = list(

    ),

    private = list(
        seqsum = NULL
    )
)
