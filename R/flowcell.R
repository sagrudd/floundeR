
#' R6 Class for performing Flowcell centric analyses
#'
#' @description
#' This class aims to simplify the handling and exploration of Flowcell based
#' data and contains various presets, designs and visualisation tools required
#' for assessing flowcell performance and metrics.
#'
#' @import R6
#'
#' @export
Flowcell <- R6::R6Class(
  inherit = FloundeR,
  classname = "Flowcell",
  public = list(

    #' @description
    #' Creates a new Flowcell object. This
    #' initialisation method performs other sanity checking of the defined
    #' file(s) and creates the required data structures.
    #'
    #' @return A new `Flowcell` object.
    initialize = function() {

    },

    #' @description
    #' set channel count summary information
    #'
    #' This method is used to provide primitive channel count information for
    #' the number of total reads that have been observed per channel - this is
    #' used for the generation of spatial plots
    #'
    #' @param channel_counts a tibble of count information
    set_channel_counts = function(channel_counts) {
      private$channel_counts <- channel_counts
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
      return(private$channel_counts)
    }
  ),


  active = list(
    #' @field platform
    #' Have a guess at the most likely flowcell platform used
    #'
    #' The sequencing summary file contains no information on the sequencing
    #' device or flowcell used. For the preparation of channel density maps
    #' it is worth considering which flowcell type is most likely to have
    #' been used - this can be guessed on the number of channels described
    #' within the data
    platform = function(value) {

      if (missing(value)) {
        if (is.null(private$platform_name)) {
          if (max(private$channel_counts$channel) < 130) {
            private$platform_name <- "Flongle"
          } else if (max(private$channel_counts$channel) > 1000) {
            private$platform_name <- "PromethION"
          } else {
            private$platform_name <- "MinION"
          }
        }
        return(private$platform_name)
      } else {
        private$platform_name <- value
      }
    },


    #' @field density_data produce channelMap for spatial plots
    #'
    #' prepares a matrix of X, Y coordinates and
    #' the corresponding readcount information for the type of flowcell
    #' predicted by `get_flowcell_platform`
    density_data = function() {
      channelMapMatrix <- reshape2::acast(
        private$.get_channel_counts(),
        col ~ row,
        value.var = "count")
      return(Angenieux$new("XYDensity", channelMapMatrix))
    }
  ),

  private = list(
    platform_name = NULL,
    channel_counts = NULL,


    .get_channel_counts = function() {
      channelMap <- private$.get_channel_map()
      channelCounts <-
        as.data.frame(matrix(rep(0, max(channelMap$channel)), ncol = 1))

      channelCounts[private$channel_counts$channel,] <- private$channel_counts$n
      channelMap <- merge(channelMap, channelCounts, by.x = "channel", by.y = 0)
      colnames(channelMap)[4] <- "count"
      return(channelMap)
    },


    .get_channel_map = function() {
      if (self$platform == "MinION") {
        return(private$.get_minion_channel_map())
      } else if (self$platform == "Flongle") {
        return(private$.get_flongle_channel_map())
      } else if (self$platform == "PromethION") {
        return(private$.get_promethion_channel_map())
      }
    },


    .get_minion_channel_map = function() {
      # build the map for R9.4.1 flowcell, as a long-form dataframe
      blockCalc <- function(i) {
        m <- matrix(seq(i, i + 63, by = 1), ncol = 8, byrow = TRUE)
        cbind(m[seq(5, 8, by = 1), ], m[seq(4), rev(seq(8))])
      }
      layout <- do.call(rbind, lapply(
        c(1, 449, 385, 321, 257, 193, 129, 65), blockCalc))
      # transpose the layout for cleaner presentation ...
      layout <- t(layout)
      channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)), which(
        layout == as.vector(layout), arr.ind = TRUE)))
      return(channelMap)
    },

    .get_flongle_channel_map = function() {
      layout <- matrix(c(seq(1, 12), 0, seq(13, 24), 0, seq(25, 114), 0,
                         seq(115, 126), 0), ncol = 13, byrow = TRUE)
      layout <- layout[rev(seq(10)), ]
      channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)),
                                        which(layout == as.vector(layout), arr.ind = TRUE)))
      return(channelMap)
    },

    .get_promethion_channel_map = function() {
      chunk <- function(i) {
        m <- matrix(seq_len(250), ncol=10, byrow=TRUE)
        m + i
      }
      layout <- do.call(cbind, lapply(seq(from=0, to=2750, by=250), chunk))
      channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)),
                                        which(layout == as.vector(layout), arr.ind = TRUE)))
      return(channelMap)
    }


  )
)
