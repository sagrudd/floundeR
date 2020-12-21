

#' R6 Class for analysing sequence sets with accompanying temporal data
#'
#' @importFrom tidyr drop_na
#' @importFrom magrittr %>%
#'
#' @export
TemporalSet <- R6::R6Class(
    inherit = FloundeR,
    classname = "TemporalSet",
    public = list(
        #' @description
        #' Initialise a new instance of the R6 Class `TemporalSet`
        #'
        #' @param seqsum a tibble of sequencing summary information
        initialize = function(seqsum=NA) {
            if (tibble::is_tibble(seqsum)) {
                if (! all(c("start_time", "duration") %in% colnames(seqsum))) {
                    stop("both [start_time] and [duration] columns required")
                }
                private$seqsum <- seqsum
                private$seqsum["speed"] <-
                    private$seqsum$sequence_length_template /
                    private$seqsum$duration
            } else {
                stop("TemporalSet requires tibble as input")
            }
        },

        #' @description
        #' Export the imported dataset(s) as a tibble
        #'
        #' @return A tibble representation of the starting dataset
        as_tibble = function() {
            return(private$seqsum)
        },

        #' @description
        #' get a rounded runtime from temporal information in sequencing summary
        #'
        #' The sequencing summary file contains information on pore dwell time
        #' and start time for given sequencing reads. An analysis of the start
        #' time information can be used to calculate a canonical sequencing
        #' run time that is rounded to key temporal breaks.
        #'
        #' @return a numeric describing the rounded runtime in hours
        runtime = function() {
            runtime <- max(private$seqsum$start_time) / 60 / 60
            expectedRuntimes <- c(4,8,12,24,36,48,64,72,96)
            temporaldistance <- sqrt((expectedRuntimes - runtime)^2)
            rruntime <- expectedRuntimes[[which(
                temporaldistance==min(temporaldistance))]]
            return(rruntime)
        },

        #' @description
        #' Prepare summary statistics that describe a flowcell run's yield
        #'
        #' A flowcell can yield both sequence bases and sequence reads and the
        #' acquisition of these data has a temporal element. This method is
        #' used to summarise run performance through assessment of yield per
        #' unit time.
        #'
        #' @param resolution describes the temporal resolution (in minutes) by
        #' which yield will be summarised.
        #' @param bases a logical that describes whether the method will
        #' summarise sequence reads or sequence bases - (TRUE) for bases by
        #' default
        #' @param cumulative defines whether the number of reads(bases) per bin
        #' is described as actual number or as a temporally cumulative number
        #'
        #' @return Angenieux object prepared for rendering in reports
        run_yield = function(resolution=15, bases=TRUE, cumulative=TRUE) {
            private$.temporal_bins(resolution)
            rls <- NULL
            if (bases) {
                rls <- private$seqsum %>%
                    dplyr::group_by(.dots=c("passes_filtering", "bin")) %>%
                    dplyr::summarize(
                        count=sum(sequence_length_template), .groups="drop")
            } else {
                rls <- private$seqsum %>%
                    dplyr::group_by(.dots=c("passes_filtering", "bin")) %>%
                    dplyr::summarize(count=dplyr::n(), .groups="drop")
            }
            rls <- tidyr::drop_na(rls)
            if (cumulative) {
                rls <- rls %>%
                    dplyr::group_by(passes_filtering) %>%
                    dplyr::arrange(bin) %>%
                    dplyr::mutate(count=cumsum(count))
            }

            Angenieux$new("2D_count", rls)
        },


        #' @description
        #' Report a binned temporal data facet (such as speed, quality, length)
        #'
        #' The temporal sequencing information can be used to summarise other
        #' data facets in a time dependent manner to address questions such as
        #' whether there is a change in sequencing speed, length or quality
        #' over time.
        #'
        #' @param resolution describes the temporal resolution (in minutes) by
        #' which yield will be summarised (60 minutes by default).
        #' @param passes is a logical that defines whether the plot should
        #' present data that has passed or failed QC - (TRUE) by default to
        #' select for only the QC passed sequence reads.
        #' @param feature defines the feature to summarise (speed by default)
        #' but could include e.g. `sequence_length_template` or
        #' `mean_qscore_template`.
        #'
        #' @return Angenieux object prepared for rendering in reports
        feature_over_time = function(
            resolution=60, passes=TRUE, feature="speed") {
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

        #' @description
        #' calculate T50 timepoint at which 50% of sequence reads are obtained
        #'
        #' Temporal data and yield data can be used to identify timepoints at
        #' which a given fraction of the data has been obtained.
        #'
        #' @param passes is a logical that defines that we should focus only
        #' on the passed QC reads (TRUE by default) - should probably be
        #' deprecated since asking for T50 of failed reads is just silly?
        #' @param t is the fractional timepoint (0.5 by default) where we are
        #' interested in knowing the time at which this fraction of reads was
        #' obtained.
        #'
        #' @return numeric describing timepoint in hours at which fractional
        #' data was obtained.
        t50 = function(passes=TRUE, t=0.5) {
            # only uses the passed data here ...
            data.filt = private$seqsum %>%
                dplyr::filter(passes_filtering == passes) %>%
                dplyr::arrange(start_time) %>%
                dplyr::summarize(
                    start_time,
                    count=cumsum(sequence_length_template), .groups="drop")
            bases = data.filt$count[length(data.filt$count)]
            message(paste0("[",bases,"] bases of sequence assessed"))
            base50.ideal <- as.numeric(bases)*t
            distance <- sqrt((as.numeric(x$count) - base50.ideal)^2)
            index <- which(distance == min(distance))
            base50 <- data.filt$count[index]
            message(
                paste0("base[",t*100,"] == ",base50," -- (",base50.ideal,")"))
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
