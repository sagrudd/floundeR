
#' @export
BamFile <- R6::R6Class(
  inherit = FloundeR,
  classname = "BamFile",
  public = list(


    initialize = function(bamfile) {
      private$bamfile <- bamfile
    }

  ),

  private = list(
    bamfile = NA

  )
)
