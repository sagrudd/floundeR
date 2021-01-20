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
    #' @param chunk_size the number of sequences to parse or process per chunk;
    #' this is used to avoid loading complete sequence collections into memory.
    #' @return A new `Fasta` object.
    #'
    #' @examples
    #' canonical_fasta <- flnDr("cluster_cons.fasta.bgz")
    #' fasta <- Fasta$new(canonical_fasta)
    initialize = function(fasta_file, chunk_size=25000) {
      private$chunk_size <- chunk_size
      cli::cli_div(theme = list(span.emph = list(color = "orange")))
      cli::cli_h1(
        stringr::str_interp(
          "creating {.emph floundR::fasta} with [${basename(fasta_file)}]"))

      if (grepl("(.fa$)|(.fasta$)", fasta_file)) {
         Rsamtools::bgzip(fasta_file, overwrite=TRUE)
        fasta_file <- sprintf("%s.bgz", fasta_file)
         cli::cli_alert(
           stringr::str_interp(
             "fasta file compressed to bzip2 format [${basename(fasta_file)}]"))
      }
      if (!any(
        grepl(
          "\\.fai$",
          list.files(dirname(fasta_file), pattern=basename(fasta_file))))) {
         cli::cli_alert_warning(
           stringr::str_interp("Creating index for [${basename(fasta_file)}]"))
         Rsamtools::indexFa(fasta_file)
      } else {
         cli::cli_alert_info(
           stringr::str_interp("index for [${basename(fasta_file)}] found"))
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
      cli::cli_alert_info(
        stringr::str_interp(
          "[${length(private$starts)}] chunks of [${chunk_size}] reads defined"))
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
        cli::cli_alert_info(
          stringr::str_interp(
            "Extracting fasta entries for sequences [${start}-${end}]"))
        dna <- Rsamtools::scanFa(
          private$fasta_file, param=private$fasta_index[start:end])
        return(dna)
      } else if (class(starts)=="numeric" & id > length(private$starts)) {
        cli::cli_alert_warning(
          stringr::str_interp(
            "Only [${length(private$starts)}] sequence chunks are available"))
        return(NA)
      } else if (is.na(private$starts)) {
        cli::cli_alert_warning(
          "please prepare sequence chunks `$sequence_chunks()` before method")
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
        dna_tib <- tibble::tibble(
          accession=names(chunk),
          width=chunk@ranges@width,
          sequence=unlist(as.character(chunk)))
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
    },

    #' @description
    #' Export the imported dataset(s) as a tibble
    #'
    #' The Fasta R6 object consumes a fasta format file and creates an object in
    #' memory that can be explored, sliced and filtered. This method dumps
    #' out the in-memory object for further exploration and development.
    #'
    #' @return A tibble representation of the starting dataset
    as_tibble = function() {
      return(private$fasta_parsed)
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
          seqsum=private$fasta_parsed)
      }
      return(private$sequencing_set)
    }

  ),

  private = list(
    sequencing_set = NULL,
    chunk_size = NULL,
    fasta_file = NULL,
    fasta_index = NULL,
    fasta_parsed = NULL,
    starts = NA,
    ends = NA,

    load_index = function() {
      cli::cli_alert_info(
        stringr::str_interp(
          "loading fasta index [${basename(private$fasta_file)}.idx]"))
      private$fasta_index <- Rsamtools::scanFaIndex(private$fasta_file)

      private$fasta_parsed <- tibble::tibble(
        read_id=as.vector(seqnames(private$fasta_index)),
        sequence_length_template=width(private$fasta_index),
        passes_filtering=TRUE)
      cli::cli_alert_success(
        stringr::str_interp(
          "[${nrow(private$fasta_parsed)}] fasta entries parsed"
        ))
    }
  )
)
