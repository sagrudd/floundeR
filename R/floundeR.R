
#' R6 Class for floundeR based analyses
#'
#' @description
#' This is an overarching class which other flounder classes are expected to
#' inherit from.
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
    }
  )
)
