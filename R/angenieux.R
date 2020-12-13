
# not aplanat (or fisheye)
#' R6 Class for visualising floundeR based datasets
#'
#' @description
#' This class aims to provide aplanat-like visualisation abstraction for the
#' floundeR framework
#'
#' @import R6
#'
#' @export
Angenieux <- R6::R6Class(
  inherit = FloundeR,
  classname = "Angenieux",
  public = list(
    
    initialize = function(key, value) {
      if (key == "XYDensity") {
        self$XYDensity(value)
      } else {
        stop(paste0("Graph type [",key,"] not implemented"))
      }
    },
    
    XYDensity = function(matrix) {
      private$graph_type = "XYDensity"
      if (!is.matrix(matrix)) {
        stop("this requires a matrix")
      }
        private$graph_data <- matrix
    }
    
  ),
  
  active = list(
    data = function(value) {
      if (missing(value)) {
        if (private$graph_type == "XYDensity") {
          return(private$graph_data)
        } else {
          stop(paste0("Graph type [",private$graph_type,"] not implemented"))
        }
      }
    }
  ),
  
  private = list(
    
    graph_type = NULL,
    graph_data = NULL
  )
)