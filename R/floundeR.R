
#' R6 Class for floundeR based analyses
#'
#' @description
#' This is an overarching class which other flounder classes are expected to
#' inherit from. The `floundeR` object contains methods and variables that are
#' maintained for other objects - there is probably not a sensible usecase for
#' actually creating a `floundeR` instance on its own.
#'
#' @import R6
#'
#' @export
FloundeR <- R6::R6Class(
  classname = "FloundeR",
  public = list(

    #' @description
    #' Creates a new FloundeR object.
    #'
    #' @return A new `FloundeR` object.
    initialize = function() {

    },

    #' @description
    #' this print method overrides the standard print function as included with
    #' R6 objects - this is to better define what is contained within an object
    #' and to provide better floundeR abstraction
    #'
    #' @param ... additional stuff passed on
    #'
    #' @return nothing (at present) - output to stdout
    print = function(...) {
      cat(paste0("<floundeR::", class(self)[1],">"))
    },


    #' @description
    #' general method for binning continuous data into sensible bins - used by
    #' various methods within the broader floundeR framework.
    #'
    #' @param data - the dataset (series) to be binned
    #' @param bins - the number of bins to prepare
    #' @param outliers - the fraction of longest outlying features that should
    #' be excluded to simplify the plot.
    #' @param qmax - a manual upper limit to use
    #'
    #' @return bin assignments for `data` elements
    bin_data = function(data, bins=20, outliers=0.025, qmax=NA) {
      if (is.na(qmax)) {
        qmax <- quantile(x=data, probs=c(1-outliers))
        qmax <- private$roundUpNice(qmax)
      }
      break_interval = private$roundUpNice(qmax / bins)
      breaks = seq(0, to = qmax, by = break_interval)
      bin_assignments <- cut(
        data, breaks, label=head(breaks, -1), include.lowest=TRUE, right=FALSE)
      return(bin_assignments)
    }
  ),

  private = list(


    # https://stackoverflow.com/questions/6461209/
    # how-to-round-up-to-the-nearest-10-or-100-or-x
    roundUpNice = function(x, nice = seq(from = 1, to = 10, by = 0.25)) {
      if (length(x) != 1)
        stop("'x' must be of length 1")
      10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
    }

  )



)
