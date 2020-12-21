
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
      },

      read_length_bins = function(qfilt=TRUE, normalised=TRUE, cumulative=FALSE, bins=20, outliers=0.025) {
          bins <- self$bin_data(
              private$seqsum$sequence_length_template, bins=bins,
              outliers=outliers)
          private$seqsum["bin"] <- as.numeric(levels(bins))[bins]

          rls <- NULL
          if (normalised) {
              rls <- private$seqsum %>%
                  dplyr::filter(.data[["passes_filtering"]]==qfilt) %>%
                  dplyr::group_by(.dots=c("barcode_arrangement", "bin")) %>%
                  dplyr::summarize(count=sum(sequence_length_template), .groups="drop")
          } else {
              rls <- private$seqsum %>%
                  dplyr::filter(.data[["passes_filtering"]]==qfilt) %>%
                  dplyr::group_by(.dots=c("barcode_arrangement", "bin")) %>%
                  dplyr::summarize(count=dplyr::n(), .groups="drop")
          }
          rls <- tidyr::drop_na(rls)

          if (cumulative) {
              rls <- rls %>% dplyr::group_by(.data[["barcode_arrangement"]]) %>% dplyr::arrange(desc(bin)) %>% dplyr::mutate(count=cumsum(count))
          }

          Angenieux$new("2D_count", rls)

      },


      quality_bins = function(bins=20, outliers=0) {
          bins <- self$bin_data(
              private$seqsum$mean_qscore_template, bins=bins,
              outliers=outliers)
          private$seqsum["bin"] <- as.numeric(levels(bins))[bins]
          rls <- private$seqsum %>%
              dplyr::group_by(.dots=c("barcode_arrangement", "bin")) %>%
              dplyr::summarize(count=dplyr::n(), .groups="drop")
          # na may creep into the data ...
          rls <- tidyr::drop_na(rls)

          Angenieux$new("2D_count", rls)
      },


      sequencingset = function(barcode) {
          return(SequencingSet$new(
              keycol = "passes_filtering",
              seqsum = seqsum %>%
                  dplyr::filter(.data[["barcode_arrangement"]]==barcode)))
      },

      temporalset = function(barcode) {
          return(TemporalSet$new(
              seqsum = seqsum %>%
                  dplyr::filter(.data[["barcode_arrangement"]]==barcode)))
      }

    ),

    active = list (

        enumerate = function(value) {
            if (missing(value)) {
                enumerated_counts <- private$seqsum %>% dplyr::group_by(.dots=c("passes_filtering", "barcode_arrangement")) %>% dplyr::summarize(count=dplyr::n(), .groups="drop")
                Angenieux$new("2D_count", enumerated_counts)
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
        private$seqsum <- dplyr::inner_join(
            seqsum, barcodedata, by="read_id") %>%
            dplyr::select(-read_id)
        return(nrow(private$seqsum) / nrow(barcodedata) * 100)
      }
    ),
)
