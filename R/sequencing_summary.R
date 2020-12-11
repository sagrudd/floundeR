
##############################################################
#' S4 sequencing summary class
#'
#' The sequencing_summary file is produced by either the Guppy or other
#' base-calling software and describes temporal and qualitative attributes of
#' the sequence reads that are produced. This class is used to process the
#' sequencing summary file and to help access results and other critical
#' information.
#'
#' @param x Description of \code{x}. The main argument in this
#'  example. Most often has such and such properties.
#'
#' @param y Description of \code{y}. An argument that is rarely
#'  used by \code{"helloworld"} methods.
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return A helloworld-ified argument. Oh, you'll see.
#'
#' @seealso \code{\link{print}} and \code{\link{cat}}
#'
#' @export
#' @docType methods
#' @rdname helloworld-methods
#'
#' @examples
#' sequencing_summary <- system.file("extdata",
#'     "sequencing_summary.txt.bz2", package="floundeR")
#' barcodes_summary <- system.file("extdata",
#'     "barcoding_summary.txt.bz2", package="floundeR")
#' seqsum <- SequencingSummary(
#'     sequ_sum=sequencing_summary,
#'     barc_sum=barcodes_summary)
SequencingSummary <- setClass(
    # Set the name for the class
    "SequencingSummary",

    # Define the slots
    slots = c(
        sequ_sum = "character",
        barc_sum = "character"
    ),

    # Set the default values for the slots. (optional)
    prototype=list(
        sequ_sum = NULL,
        barc_sum = NULL
    ),
    # Make a function that can test to see if the data is consistent.
    # This is not called if you have an initialize function defined!
    validity=function(object)
    {
        if (!file.exists(object@sequ_sum)) {
            return(paste0("sequ_sum [",path,"] does not exist"))
        }

        if (!is.null(object@barc_sum) && !file.exists(object@barc_sum)) {
            return(paste0("barc_sum [",path,"] does not exist"))
        }
        return(TRUE)
    }
)






setGeneric("helloworld", function(x, y, ...){
    cat("Hello World!")
    cat("\n")
    standardGeneric("helloworld")
})

#' @rdname helloworld-methods
#' @aliases helloworld,ANY,ANY-method
setMethod("helloworld", "ANY", function(x, y, ...){
    cat(class(x))
})

#' @rdname helloworld-methods
#' @aliases helloworld,character,ANY-method
setMethod("helloworld", "character", function(x){
    show(x)
})

#' @rdname helloworld-methods
#' @aliases helloworld,character,character-method
setMethod("helloworld", c("character", "character"), function(x, y){
    show(x)
})
