

#' R6 Class for loading and analysing nanopore (and other) FASTQ files
#'
#' @description
#' This class aims to simplify the handling and exploration of
#' FASTQ files and provides simple methods for accessing
#' information that can be used to assess the contents of a FASTQ file.
#'
#' @import R6
#' @importFrom magrittr %>%
#'
#' @export
Fastq <- R6::R6Class(
  inherit = FloundeR,
  classname = "Fastq",
  public = list(
    #' @field fastq_file the file.path
    #' to the query FAST5 file
    fastq_file = NULL,
    
    #' @description
    #' Creates a new SequencingSummary object. This
    #' initialisation method performs other sanity checking
    #' of the defined file(s) to ensure that it is indeed
    #' parseable and creates the required data structures.
    #'
    #' @param fastq_file The source
    #' sequencing_summary file.
    #' @return A new `Fastq` object.
    #'
    #' @examples
    #' canonical_fastq <- flnDr("example.fastq.gz")
    #' fastq <- Fastq$new(canonical_fastq)
    initialize = function(
      fastq_file) {
      
      # first pass QC - let's ensure that this is a single FILE
      .check_path(
        "fastq_file", fastq_file)
      self$fastq_file = fastq_file
      
      # and import the fastq metadata (cleanly)
      if (!private$.parse_fastq()) {
        stop("Failed to import the fastq file provided ...")
      }
    },
    
    #' @description
    #' Export the imported dataset(s) as a tibble
    #'
    #' This object consumes a sequencing summary file (and optionally the
    #' corresponding barcoding_summary file) and creates an object in
    #' memory that can be explored, sliced and filtered. This method dumps
    #' out the in-memory object for further exploration and development.
    #'
    #' @return A tibble representation of the starting dataset
    as_tibble = function() {
      return(private$fastq)
    }
    
  ),
  
  active = list(
    
    #' @field sequencingset
    #' The `sequencingset` active binding returns a sequencingset object
    #' that is canonically structured around the `passes_filtering` logical
    #' field to allow assessment of sequencing characteristics.
    sequencingset = function(column="passes_filtering", keys=NULL) {
      if (is.null(private$sequencing_set)) {
        private$sequencing_set <- SequencingSet$new(
          keycol=column,
          seqsum=private$fastq)
      }
      return(private$sequencing_set)
    }
    
  ),
  
  
  private = list(
      fastq = NULL,
      sequencing_set = NULL,
      
      .parse_fastq = function() {
        private$fastq <- tibble::as_tibble(fishy_fastq(self$fastq_file))
        return(TRUE)
      }
    
    
  )
)
