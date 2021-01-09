# canonical_fasta <- flnDr("cluster_cons.fasta.bgz")
# fasta <- Fasta$new(canonical_fasta)
# x <- fasta$sequence_chunks()
# fasta$get_sequence_chunk(1)

#' R6 Class for loading and analysing FASTA files
#'
#' @description
#' This class aims to simplify the handling and exploration of FASTA files and
#' provides simple methods for accessing information that can be used to assess
#' bulk contents from a FASTA file - the analysis framework is provided by
#' Rsamtools (alone).
#'
#' @import R6
#' @importFrom magrittr %>%
#'
#' @export
Fasta <- R6::R6Class(
  inherit = FloundeR,
  classname = "Fasta",
  public = list(

    #' @description
    #' Creates a new Fasta object. This
    #' initialisation method performs other sanity checking
    #' of the defined file(s) to ensure that it is indeed
    #' parseable and creates the required data structures.
    #'
    #' @param fasta_file The source
    #' sequencing_summary file.
    #' @return A new `Fastq` object.
    #'
    #' @examples
    #' canonical_fasta <- flnDr("cluster_cons.fasta.bgz")
    #' fasta <- Fasta$new(canonical_fasta)
    initialize = function(fasta_file) {
      cli::cli_div(theme = list(span.emph = list(color = "orange")))
      cli::cli_h1(stringr::str_interp("creating {.emph floundR::fasta} with [${basename(fasta_file)}]"))

      if (grepl("(.fa$)|(.fasta$)", fasta_file)) {
         Rsamtools::bgzip(fasta_file, overwrite=TRUE)
        fasta_file <- sprintf("%s.bgz", fasta_file)
         cli::cli_alert(stringr::str_interp("fasta file compressed to bzip2 format [${basename(fasta_file)}]"))
      }
      if (!any(grepl("\\.fai$", list.files(dirname(fasta_file), pattern=basename(fasta_file))))) {
         cli::cli_alert_warning(stringr::str_interp("Creating index for [${basename(fasta_file)}]"))
         Rsamtools::indexFa(fasta_file)
      } else {
         cli::cli_alert_info(stringr::str_interp("index for [${basename(fasta_file)}] found"))
      }
      private$fasta_file <- fasta_file

      private$load_index()
    },

    #' @description
    #' Split the fasta sequence file explored by the package into sequence
    #' chunks for e.g. import into a relational database.
    #' 
    #' @param chunk_size The number of fasta entries that should be contained
    #' within a single chunk (default: 1000)
    #' 
    #' @return an invisible integer that defines the number of possible chunks;
    #' this can for example be iterated over
    sequence_chunks = function(chunk_size=1000) {
      max <- length(private$fasta_index)
      private$starts <- seq(from=1, to = max, by = chunk_size)
      private$ends <- c(seq(from=0, to = max, by = chunk_size)[-1], max)
      cli::cli_alert_info(stringr::str_interp("sequence collection can be resolved into [${length(private$starts)}] chunks of [${chunk_size}] reads"))
      invisible(length(private$starts))
    },

    #' @description 
    #' Get a chunk of fasta sequences from a larger monolithic file
    #' 
    #' @param id the chunk (see `$sequence_chunks()`) to extract sequence for -
    #' this must be an integer that is > 0 and <= sequence_chunks.
    #' 
    #' @return DNAStringSet containing the fasta entries corresponding to the
    #' specified sequence chunk.
    get_sequence_chunk = function(id=1) {
      if (class(private$starts)=="numeric" & id <= length(private$starts)) {
        start = private$starts[id]
        end = private$ends[id]
        cli::cli_alert_info(stringr::str_interp("Extracting fasta entries for sequences [${start}-${end}]"))
        dna <- Rsamtools::scanFa(private$fasta_file, param=private$fasta_index[start:end])
        return(dna)
      } else if (class(starts)=="numeric" & id > length(private$starts)) {
        cli::cli_alert_warning(stringr::str_interp("Only [${length(private$starts)}] sequence chunks are available"))
        return(NA)
      } else if (is.na(private$starts)) {
        cli::cli_alert_warning("please prepare sequence chunks `$sequence_chunks()` before this method")
        return(NA)
      }
    },

    #' @description 
    #' Get a chunk of fasta sequences from a larger monolithic file as a tibble
    #'
    #' @param id the chunk (see `$sequence_chunks()`) to extract sequence for -
    #' this must be an integer that is > 0 and <= sequence_chunks.
    #' 
    #' @return `tibble` containing the fasta entries corresponding to the
    #' specified sequence chunk.
    get_tibble_chunk = function(id) {
      chunk <- self$get_sequence_chunk(id)
      if ("DNAStringSet" %in% class(chunk)) {
        dna_tib <- tibble::tibble(accession=names(chunk), width=chunk@ranges@width, sequence=unlist(as.character(chunk)))
        return(dna_tib)
      }
      return(chunk)
    },
    
    #' @description 
    #' return the Rsamtools FASTA index
    #' @return GRanges object describing the fasta elements contained within the
    #' sequence file
    get_index = function() {
      return(private$fasta_index)
    },
    
    #' @description 
    #' return the number of sequence elements contained within the sequence
    #' file specified
    #' @return integer of fasta entries in file
    count = function() {
      return(length(private$fasta_index))
    }

  ),

  private = list(
    fasta_file = NULL,
    fasta_index = NULL,
    starts = NA,
    ends = NA,

    load_index = function() {
      cli::cli_alert_info(stringr::str_interp("loading fasta index [${basename(private$fasta_file)}.idx]"))
      private$fasta_index <- Rsamtools::scanFaIndex(private$fasta_file)
    }
  )
)
