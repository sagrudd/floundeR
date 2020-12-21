
#' R6 Class for loading, visualising and analysing barcode information
#'
#' @description
#'
#'
#' @export
MultiplexSet <- R6::R6Class(
    inherit = FloundeR,
    classname = "MultiplexSet",
    public = list(
      initialize = function(seqsum=NA, barcoding_summary_file=NA) {
        if (tibble::is_tibble(seqsum)) {
          if (!"barcode_arrangement" %in% colnames(seqsum) &
              (is.na(barcoding_summary_file) || is.null(barcoding_summary_file))) {
            stop("MultiplexSet requires a `barcode_arrangement` column")
          } else {
            imported <- private$.import_barcoding_summary(
                seqsum, barcoding_summary_file)
            print(imported)
            if (imported < 100) {
            stop(paste0(
                "MultiplexSet failed trying to import barcoding file - ",
                "[",imported,"]% of reads matched seq/barc files"))
            }
          }
        } else {
          stop("MultiplexSet requires tibble as constructor")
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
          return(private$seqsum)
      }

    ),

    private = list(
      seqsum = NULL,

      .import_barcoding_summary = function(seqsum, barcoding_summary_file) {
        barcodedata <- readr::read_tsv(
          barcoding_summary_file,
          col_types = readr::cols_only(
            "read_id"="c",
            "barcode_arrangement"="c"))
        private$seqsum <- dplyr::inner_join(
            seqsum, barcodedata, by="read_id") %>%
            dplyr::select(-read_id)
        return(nrow(private$seqsum) / nrow(barcodedata) * 100)
      }
    ),
)
