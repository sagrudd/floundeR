
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
                private$seqsum["speed"] <-
                    private$seqsum$sequence_length_template / private$seqsum$duration
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

        run_yield = function(resolution=15, bases=TRUE, cumulative=TRUE) {
            private$.temporal_bins(resolution)
            rls <- NULL
            if (bases) {
                rls <- private$seqsum %>%
                    dplyr::group_by(.dots=c("passes_filtering", "bin")) %>%
                    dplyr::summarize(count=sum(sequence_length_template), .groups="drop")
            } else {
                rls <- private$seqsum %>%
                    dplyr::group_by(.dots=c("passes_filtering", "bin")) %>%
                    dplyr::summarize(count=dplyr::n(), .groups="drop")
            }
            rls <- tidyr::drop_na(rls)
            if (cumulative) {
                rls <- rls %>% dplyr::group_by(passes_filtering) %>% dplyr::arrange(bin) %>% dplyr::mutate(count=cumsum(count))
            }

            Angenieux$new("2D_count", rls)
        },


        feature_over_time = function(resolution=60, passes=TRUE, feature="speed") {
            private$.temporal_bins(resolution)

            if (!feature %in% colnames(private$seqsum)) {
                stop("feature_over_time requires a valid seqsum column name")
            }
            data.filt = private$seqsum %>% filter(passes_filtering == x)
            data.tibble = tibble::tibble(
                data.filt$bin,
                data.filt[feature])
            colnames(data.tibble) <- c("bin", feature)
            Angenieux$new("boxplot", data.tibble)
        },


        t50 = function(passes=TRUE, t=0.5) {
            # only uses the passed data here ...
            data.filt = private$seqsum %>%
                dplyr::filter(passes_filtering == passes) %>%
                dplyr::arrange(start_time) %>%
                dplyr::summarize(start_time, count=cumsum(sequence_length_template), .groups="drop")
            bases = data.filt$count[length(data.filt$count)]
            message(paste0("[",bases,"] bases of sequence assessed"))
            base50.ideal <- as.numeric(bases)*t
            distance <- sqrt((as.numeric(x$count) - base50.ideal)^2)
            index <- which(distance == min(distance))
            base50 <- data.filt$count[index]
            message(paste0("base[",t*100,"] == ",base50," -- (",base50.ideal,")"))
            return(data.filt$start_time[index] / 60 / 60)

        }


    ),

    active = list(

    ),

    private = list(
        seqsum = NULL,
        .temp_resolution = -1,

        .temporal_bins = function(resolution) {
            if (resolution == private$.temp_resolution) {
                return()
            }
            bins <- self$bin_data(
                private$seqsum$start_time/60/60,
                bins=self$runtime()*60/resolution,
                qmax=self$runtime())
            # bins is thus in hours with resolution in minutes ...
            private$seqsum["bin"] <- as.numeric(levels(bins))[bins]
            private$.temp_resolution <- resolution
        }
    )
)
