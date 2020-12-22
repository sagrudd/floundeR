
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
#' @importFrom tools file_ext
#'
#' @export
Angenieux <- R6::R6Class(
  inherit = FloundeR,
  classname = "Angenieux",
  public = list(

    #' @field colourMap - default colourMap for plots requiring discrete colours
    colourMap = NA,

    #' @description
    #' Creates a new Angenieux object. This
    #' initialisation method performs other sanity checking of the defined
    #' file(s) and creates the required data structures.
    #'
    #' @param key the datatype that is being passed to the object
    #'
    #' @param value the data that is being passed to the object
    #'
    #' @return A new `Angenieux` object.
    initialize = function(key, value) {

      self$colourMap = RColorBrewer::brewer.pal(8, "Set1")

      private$hm.palette = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Blues"), space = "Lab")
      private$.plot_elements <- list()

      if (key == "XYDensity") {
        if (!is.matrix(value)) {
          stop("this requires a matrix")
        }
        private$graph_type = key
        private$graph_data <- value
      } else if (key == "1D_count") {
        if (!tibble::is_tibble(value)) {
          stop("this requires a tibble")
        }
        private$graph_type = key
        private$graph_data <- value
      } else if (key == "2D_count") {
        if (!tibble::is_tibble(value)) {
          stop("this requires a tibble")
        }
        private$graph_type = key
        private$graph_data <- value
      }  else if (key == "boxplot") {
        if (!tibble::is_tibble(value)) {
          stop("this requires a tibble")
        }
        private$graph_type = key
        private$graph_data <- value
      }

      else {
        stop(paste0("Graph type [",key,"] not implemented"))
      }
    },


    #' @description
    #' add an `AngenieuxDecoration` to the plot.
    #'
    #' @param item an `AngenieuxDecoration`
    add = function(item) {
      if (!class(item)[1] == "AngenieuxDecoration") {
        stop("Can only add [AngenieuxDecoration] elements")
      }
      private$.plot_elements <- append(private$.plot_elements, item)
      invisible(self)
    },


    #' @description
    #' Prepare and present an Angenieux plot
    #'
    #' @param ... parameters passed on to downstream methods - please see
    #' examples for further examples as to how Angenieux plots can be customised
    #' using this approach.
    plot = function(...) {
      if (private$graph_type == "XYDensity") {
        return(private$.plot_xy_density(...))
      } else if (private$graph_type == "1D_count") {
        return(private$.plot_1d_count(...))
      } else if (private$graph_type == "2D_count") {
        return(private$.plot_2d_count(...))
      } else if (private$graph_type == "boxplot") {
        return(private$.plot_boxplot(...))
      } else {
        stop(paste0("Graph type [",private$graph_type,"] not implemented"))
      }
    },

    #' @description
    #' Specify that Angenieux plot should be saved to file
    #'
    #' When working at the console an Angenieux plot may be plotted directly
    #' to the console. When preparing reports through Rmarkdown or Pkgdown a
    #' more logical saving of plots to a discrete file location may make more
    #' sense. The method is used to instruct Angenieux that the plot should be
    #' saved to a given location and with a given file format.
    #'
    #' @param target_file the file with extension e.g. `figure1.png`
    #' @param width the width of figure to save (12 by default)
    #' @param height the height of figure to save (7.5 by default)
    #' @param units the unit to use for height and width (cm by default)
    #' @param dpi the plot resolution (print/300 by default)
    #'
    #' @return the original Angenieux object (self)
    to_file = function(target_file, width=18, height=12, units="cm", dpi="print") {
      private$target_type <- tolower(tools::file_ext(target_file))
      private$target_file <- target_file
      private$target_file_width = width
      private$target_file_height = height
      private$target_file_units = units
      private$target_file_dpi = dpi
      invisible(self)
    },

    #' @description
    #' Set the title used in the given Angenieux plot
    #'
    #' @param title - the title to use on the plot
    #'
    #' @return the original Angenieux object (self)
    set_title = function(title) {
      private$graph_title <- title
      invisible(self)
    }

  ),

  active = list(
    #' @field data
    #' A method to dump out the stored data from an `Angenieux` object
    data = function(value) {
      if (missing(value)) {
        if (private$graph_type == "XYDensity") {
          return(private$graph_data)
        } else if (private$graph_type == "1D_count") {
          return(private$graph_data)
        } else if (private$graph_type == "2D_count") {
          return(private$graph_data)
        } else if (private$graph_type == "boxplot") {
          return(private$graph_data)
        } else {
          stop(paste0("Graph type [",private$graph_type,"] not implemented"))
        }
      }
    }
  ),

  private = list(
    .plot_elements = NA,
    graph_type = NULL,
    graph_data = NULL,
    graph_title = "angenieux plot",
    hm.palette = NULL,
    target_type = NA,
    target_file = NA,
    target_file_width = NA,
    target_file_height = NA,
    target_file_units = NA,
    target_file_dpi = NA,

    .plot_xy_density = function(count=FALSE) {

      molten_matrix <- reshape2::melt(
        private$graph_data, value.name = "Count", varnames=c('X', 'Y'))

      plot <- ggplot2::ggplot(
        molten_matrix,
        ggplot2::aes_string(x = "X", y = "Y", fill = "Count")) +
        ggplot2::geom_tile() +
        ggplot2::scale_x_discrete(breaks = NULL) +
        ggplot2::scale_y_discrete(breaks = NULL) +
        ggplot2::coord_equal() +
        ggplot2::scale_fill_gradientn(colours = private$hm.palette(100)) +
        ggplot2::scale_color_gradient2(low = private$hm.palette(100)[100], high = private$hm.palette(100)[1]) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::theme(panel.border = element_blank(), panel.grid.major=element_blank(),
              panel.grid.minor = element_blank(), axis.title.x = element_blank(),
              axis.title.y = element_blank(), legend.position = "bottom")
      if (count) {
        plot <- plot + ggplot2::geom_text(
          data = molten_matrix,
          ggplot2::aes_string(x = "X", y = "Y", label = "Count", color = "Count"),
          show.legend = FALSE, size = 2.5)
      }
      return(private$.decorate_plot(plot))
    },


    .plot_1d_count = function(style="histogram") {
      key <- colnames(private$graph_data)[[1]]
      if (style == "stacked") {
        plot <- ggplot2::ggplot(
          private$graph_data,
          aes_string(x=as.factor(" "), y="count", fill=key)) +
          ggplot2::geom_col(width=0.2) +
          ggplot2::coord_flip()
        return(private$.decorate_plot(plot))
      } else {
        plot <- ggplot2::ggplot(
          private$graph_data,
          ggplot2::aes_string(key, "count")) +
          ggplot2::geom_bar(stat = "identity", width = 0.5)
        return(private$.decorate_plot(plot))
      }
    },

    .plot_2d_count = function(style="line") {
      key <- colnames(private$graph_data)[[1]]
      level <- colnames(private$graph_data)[[2]]
      molten <- reshape2::melt(
        private$graph_data, id.vars=c(level, key), measure.vars=c("count"))
      if (style == "line") {
        plot <- ggplot2::ggplot(
          molten,
          ggplot2::aes_string(x=level, y="value", colour=key)) +
          ggplot2::geom_line()
        return(private$.decorate_plot(plot))
      } else {
        #molten[[level]] <- factor(
        #  molten[[level]],
        #  sort(unique(molten[[level]])))
        plot <- ggplot2::ggplot(
          molten,
          ggplot2::aes_string(x=level, y="value", fill=key)) +
          ggplot2::geom_bar(stat="identity") +
          ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))
        return(private$.decorate_plot(plot))
      }
    },

    .plot_boxplot = function() {
      key <- colnames(private$graph_data)[[2]]
      plot <- ggplot2::ggplot(
        private$graph_data,
        ggplot2::aes_string(x="bin", y=key, group="bin")) +
        ggplot2::geom_boxplot(fill="steelblue", outlier.shape=NA) +
        ggplot2::scale_x_continuous()
      return(private$.decorate_plot(plot))
    },


    .decorate_plot = function(plot) {
      for (decoration in private$.plot_elements) {
        plot <- plot + decoration$decoration
      }
      plot <- plot +
        ggplot2::labs(title = private$graph_title) +
        ggplot2::theme(text = element_text(size = 10))
      return(private$.handle_plot_logistsics(plot))
    },


    .handle_plot_logistsics = function(plot) {
      if (is.na(private$target_type)) {
        return(plot)
      } else if (private$target_type == "png") {
        message("saving plot as [png]")
        ggplot2::ggsave(
          private$target_file,
          plot = plot,
          device = private$target_type,
          width = private$target_file_width,
          height = private$target_file_height,
          units = private$target_file_units,
          dpi = private$target_file_dpi)
        return(private$target_file)
      } else {
        stop(paste0("plottype - [",private$target_type,"] is not defined"))
      }
    }
  )
)





#' R6 Class for describing additional Angenieux decorations
#'
#' @export
AngenieuxDecoration <- R6::R6Class(
  inherit = FloundeR,
  classname = "AngenieuxDecoration",
  public = list(

    #' @field decoration
    #' This field contains the decoration that will be applied to the Angenieux
    #' object
    decoration = NA,

    #' @description
    #' This is the constructor for the AngenieuxDecoration object
    #'
    #' @param decoration_type
    #' This field is used to specify the type of decoration; the cleanest type
    #' at the moment is currently the `ggplot2` type.
    #' @param ... the other variables passed on to methods contained within the
    #' object
    initialize = function(decoration_type, ...) {
      if (decoration_type=="vline") {
        self$.add_vline(...)
      } else if (decoration_type=="vlegend") {
        self$.add_vlegend(...)
      } else if (decoration_type=="ggplot2") {
        self$.add_ggplot2(...)
      }
    },

    #' @description
    #' Add a vertical line to a ggplot2 graph within Angenieux
    #'
    #' @param xintercept the point at which the vertical line will intercept the
    #' x-axis
    #' @param colour the colour of the line
    #' @param size the width of the line (default 1)
    .add_vline = function(xintercept, colour="green", size=1) {
      self$decoration <- ggplot2::geom_vline(
        xintercept=xintercept,
        colour=colour,
        size=size)
    },

    #' @description
    #' Add a legend text to accompany a vertical line
    #'
    #' @param xintercept the point at which the vertical line will intercept the
    #' x-axis
    #' @param colour the colour of the line
    #' @param legend the text to display at the given location
    #' @param hjust horizonal justify (0=left, 1=right)
    #' @param vjust vertical justify (0=bottom, 1=top)
    #' @param size the size of the font to present at the given location
    .add_vlegend = function(xintercept, colour="green", legend="", hjust=0, vjust=1, size = 6) {
      self$decoration <- ggplot2::annotate(
        "text",
        x=xintercept,
        y=+Inf,
        label=legend,
        hjust=hjust,
        vjust=vjust,
        colour=colour,
        size=size)
    },

    #' @description
    #' Just add some plain `ggplot2` to the AngenieuxDecoration and layer on to
    #' the Angenieux plot - this is for the lazy hacking out and visualisation
    #' of plots
    #'
    #' @param facet the stuff to be layered onto the graph
    .add_ggplot2 = function(facet) {
      self$decoration <- facet
    }
  )
)
