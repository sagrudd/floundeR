# blast_results <- flnDr("drosophila_uniref100.blastx.gz")
# blast <- Blast$new(blast_file=blast_results)

#' R6 Class for loading and analysing BLAST results in basic `Pairwise` format
#'
#' @description
#' BLAST results are a fundamental unit of basic comparative genomics. This
#' R6 object has been implemented for the systematic exploration of BLAST
#' results within the scope of the Lodestar project.
#'
#' @import R6
#' @importFrom magrittr %>%
#'
#' @export
Blast <- R6::R6Class(
  inherit = FloundeR,
  classname = "Blast",
  public = list(
    #' @field blast_file the file.path to the query BLAST results file
    blast_file = NA,

    #' @description
    #' Creates a new Blast object. This initialisation method performs other
    #' sanity checking of the defined file(s) to ensure that it is indeed
    #' parseable and creates the required data structures.
    #'
    #' @param blast_file The source
    #' sequencing_summary file.
    #' @return A new `Blast` object.
    #'
    #' @examples
    #' blast_results <- flnDr("drosophila_uniref100.blastx.gz")
    #' blast <- Blast$new(blast_file=blast_results)
    initialize = function(blast_file) {
      # first pass QC - let's ensure that this is a single FILE
      .check_path(
        "blast_file", blast_file)
      self$blast_file = blast_file
      private$parse_blast_results()
    },


    #' @description
    #' Return the number of BLAST results that are contained within the BLAST
    #' file provided.
    count = function() {
      return(private$blast_count)
    }
  ),

  private = list(

    blast_count = NA,

    parse_blast_results = function() {
      cli::cli_alert(stringr::str_interp("Parsing BLAST result file [${basename(self$blast_file)}]"))
    }
  )
)

