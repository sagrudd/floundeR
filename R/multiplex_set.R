
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
            private$.import_barcoding_summary(seqsum, barcoding_summary_file)
            stop("MultiplexSet failed trying to import barcoding file")
          }
        } else {
          stop("MultiplexSet requires tibble as constructor")
        }
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
        
        pso <- order(seqsum$read_id, method = "radix")
        seqsum <- seqsum[pso, ]
        
        bco <- order(barcodedata$read_id, method = "radix")
        barcodedata <- barcodedata[bco, ]
        
        barcodeMapping <- fmatch(seqsum$read_id, barcodedata$read_id)
        seqsum$barcode_arrangement <- barcodedata[barcodeMapping,
                                                  c("barcode_arrangement")]        
      }
    ),
)
