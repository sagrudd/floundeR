

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
    #' to the query FASTQ file
    fastq_file = NULL,

    #' @description
    #' Creates a new Fastq object. This
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
      if (!private$load_index()) {
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
      return(private$fastq_index)
    },
    
    #' @description
    #' Split the fastq sequence file explored by the package into sequence
    #' chunks for e.g. import into a relational database.
    #' 
    #' @param chunk_size The number of fastq entries that should be contained
    #' within a single chunk (default: 10000)
    #' 
    #' @return an invisible integer that defines the number of possible chunks;
    #' this can for example be iterated over
    sequence_chunks = function(chunk_size=10000) {
      private$chunk_size <- chunk_size
      max_val <- base::nrow(private$fastq_index)
      # print(paste0("max=", max))
      private$starts <- seq(from=1, to = max_val, by = chunk_size)
      # print(starts)
      private$ends <- c(seq(from=0, to = max_val, by = chunk_size)[-1], max)
      cli::cli_alert_info(
        stringr::str_interp(
          "sequence collection can be resolved into [{length(private$starts)}] chunks of [{chunk_size}] reads"))
      invisible(length(private$starts))
    },
    
    
    #' @description 
    #' Get a chunk of fastq sequences from a larger monolithic file. This method
    #' can be called for up to `$sequence_chunks()` times or until NULL results
    #' are returned.
    #' 
    #' @return tibble containing the fastq entries corresponding to the
    #' available sequence chunk.
    get_sequence_chunk = function() {
      if (is.null(private$chunk_size)) {
        cli::cli_alert_danger("Please define chunking with `$sequence_chunks()` first")
        silent_stop()
      }
      if (is.null(private$fqstream)) {
        cli::cli_alert("opening fastq stream")
        private$fqstream <- ShortRead::FastqStreamer(self$fastq_file, private$chunk_size)
      }
      
      fq <- ShortRead::yield(private$fqstream)
      if (length(fq) == 0) {
        cli::cli_alert("end of fastq file reached; closing connection")
        close(private$fqstream)
        private$fqstream <- NULL
        return(NULL)
      }
      quality <- Biostrings::quality(fq)
      qscore <- unlist(
        lapply(as.character(quality@quality), qualToMeanQ))
      fq_tib <- tibble::tibble(
        accession=as.character(ShortRead::id(fq)), 
        sequence=as.character(ShortRead::sread(fq)),
        sequence_length_template=Biostrings::width(fq), 
        quality=as.character(quality@quality),
        mean_qscore_template=qscore, 
        passes_filtering=TRUE)
      return(fq_tib)
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
          seqsum=private$fastq_index)
      }
      return(private$sequencing_set)
    }

  ),


  private = list(
      fastq_index = NULL,
      sequencing_set = NULL,
      fqstream = NULL,
      starts = NULL,
      ends = NULL,
      chunk_size = NULL,

      load_index = function(n=10000) {
        #private$fastq <- tibble::as_tibble(fishy_fastq(self$fastq_file))
        
        # DataFrame df = DataFrame::create( 
        #   Named("id") = id_vect , 
        #   _["sequence_length_template"] = length_vect ,
        #   _["mean_qscore_template"] = quality_vect ,
        #   _["passes_filtering"] = qc_pass_vect);
        if (is.null(private$fqstream)) {
          cli::cli_alert("opening fastq stream")
          private$fqstream <- ShortRead::FastqStreamer(self$fastq_file, n)
          repeat {
            fq <- ShortRead::yield(private$fqstream)
            if (length(fq) == 0) {
              break
            }
            quality <- Biostrings::quality(fq)
            qscore <- unlist(
              lapply(as.character(quality@quality), qualToMeanQ))
            fq_tib <- tibble::tibble(
              sequence_length_template=Biostrings::width(fq), 
              mean_qscore_template=qscore, 
              passes_filtering=TRUE)
            private$fastq_index <- private$fastq_index %>% dplyr::bind_rows(fq_tib)
          }
        }
        close(private$fqstream)
        private$fqstream <- NULL
        return(TRUE)
      }


  )
)
