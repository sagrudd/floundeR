
# not aplanat (or fisheye)
#' R6 Class for visualising floundeR based datasets
#'
#' @description
#' This class aims to provide aplanat-like visualisation abstraction for the
#' floundeR framework
#'
#' @import R6
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
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
    },
    
    
    plot = function(value) {
      if (missing(value)) {
        if (private$graph_type == "XYDensity") {
          return(private$.plot_xy_density())
        } else {
          stop(paste0("Graph type [",private$graph_type,"] not implemented"))
        }
      }
    }
  ),
  
  private = list(
    
    graph_type = NULL,
    graph_data = NULL,
    graph_title = "angenieux plot",
    
    .plot_xy_density = function() {
      
      hm.palette <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Blues"), space = "Lab")
      
      molten_matrix <- reshape2::melt(
        private$graph_data, value.name = "Count", varnames=c('X', 'Y'))
      
      plot <- ggplot2::ggplot(
        molten_matrix, 
        ggplot2::aes_string(x = "X", y = "Y", fill = "Count")
      ) + ggplot2::geom_tile() + ggplot2::scale_x_discrete(breaks = NULL) + 
        ggplot2::scale_y_discrete(breaks = NULL) + 
        ggplot2::coord_equal() + 
        ggplot2::scale_fill_gradientn(colours = hm.palette(100)) + 
        ggplot2::scale_color_gradient2(low = hm.palette(100), high = hm.palette(1)) + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) + 
        ggplot2::labs(title = private$graph_title) +
        ggplot2::theme(panel.border = element_blank(), panel.grid.major=element_blank(),
              panel.grid.minor = element_blank(), axis.title.x = element_blank(),
              axis.title.y = element_blank(), legend.position = "bottom",
              legend.key.width = unit(5.6, "cm")) +
        ggplot2::geom_text(
          data = molten_matrix, 
          ggplot2::aes_string(x = "X", y = "Y", label = "Count", color = "Count"),
          show.legend = FALSE, size = 2.5)
      return(plot)
    }
  )
)