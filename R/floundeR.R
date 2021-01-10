
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
    },


    #' @description
    #' scale numerics in bases into kilo/mega/giga etc measurements for e.g.
    #' plots and other graphic visualisations
    #'
    #' @param val a numeric (raw)
    #'
    #' @return character representation scaled accordingly
    num_scale = function(val) {
      if (is.na(val)) {
        return(val)
      }
      val = as.numeric(val)
      if (val > 1e+12) {
        return(paste0(val/1e+12,"Tb"))
      } else if (val > 1e+9) {
        return(paste0(val/1e+9,"Gb"))
      } else if (val > 1e+6) {
        return(paste0(val/1e+6,"Mb"))
      } else if (val > 1e+3) {
        return(paste0(val/1e+3,"kb"))
      } else {
        return(val)
      }
    },

    #' @description
    #' Apply the `num_scale` method across a multi-member vector of numerics
    #'
    #' @param vals the vector of numerics
    nums_scale = function(vals) {
      unlist(lapply(vals, self$num_scale))
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



#' calculate mean Phred scores from list of Q values
#'
#' This method calculates the mean quality value from a list of q-values;
#' handles the challenge of relative log-scaling of the Q and challenges of a
#' more linear mean; this avoids an artificial inflation of phred-scores
#'
#' @param q is a vector of q values
#' @return mean phred scaled q-value
#'
#' @examples
#' mean(c(20,30,40))
#' phredmean(c(20,30,40))
#'
#' @export
phredmean <- function(q) {
  -10 * log10(mean(10^(q/-10), na.rm = TRUE))
}



#' calculate mean Phred score from an ASCII encoded phred string
#'
#' FASTQ and BAM store per base qualities as an ASCII string. Accessory methods
#' in e.g. ShortRead allow for a sum of the numeric encoded scores; this is not
#' corrected for the log/linear so scores are synthetically boosted
#' - this simple method performs mean on the character level data ...
#'
#' @param qstr is an ASCII encoded Phred quality score
#' @return mean phred scaled q-value
#'
#' @examples
#' qualToMeanQ('ABCDEF')
#'
#' @export
qualToMeanQ <- function(qstr) {
  baseq <- as.numeric(charToRaw(qstr)) - 33
  meanerror <- mean(10^(baseq/-10))
  -10 * log10(meanerror)
}

silent_stop <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop("\r ", call.=FALSE)
}
