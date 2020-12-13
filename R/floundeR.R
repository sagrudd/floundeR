
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
    initialize = function() {
      
    },
    
    print = function(...) {
      cat(paste0("<floundeR::", class(self)[1],">"))
    }
  )
)