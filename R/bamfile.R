
#' R6 class for legacy BAM file references
#'
#' @description
#' `BamFile` is a minimal legacy holder for a BAM file path. The rebooted BAM
#' QC surface is expected to move to curated in-process Rust bindings through
#' Bamana, but this class remains exported until that replacement API is
#' defined.
#'
#' @export
BamFile <- R6::R6Class(
  inherit = FloundeR,
  classname = "BamFile",
  public = list(


    #' @description
    #' Create a new `BamFile` object.
    #'
    #' @param bamfile Path to a BAM file.
    #' @return A new `BamFile` object.
    initialize = function(bamfile) {
      private$bamfile <- bamfile
    }

  ),

  private = list(
    bamfile = NA

  )
)
