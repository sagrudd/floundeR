
#' R6 Class for loading, visualising and analysing barcode information
#'
#' @importFrom magrittr %>%
#'
#'
#' @export
MultiplexSet <- R6::R6Class(
  inherit = FloundeR,
  classname = "MultiplexSet",
  public = list(

    #' @description
    #' Initialise a new instance of the R6 Class `MultiplexSet`
    #'
    #' @param seqsum a tibble of sequencing summary information
    #' @param barcoding_summary_file is a file.path to the corresponding
    #' barcoding_summary file that should be merged with the `seqsum`
    #' content.
    initialize = function(seqsum=NA, barcoding_summary_file=NA) {
      if (tibble::is_tibble(seqsum)) {
        if (!"barcode_arrangement" %in% colnames(seqsum) &
            (is.na(barcoding_summary_file) ||
             is.null(barcoding_summary_file))) {
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


    #' @description
    #' bin the sequences in `seqsum` content into bins of sequence length
    #'
    #' The nanopore sequencing run is expected to return a collection of
    #' sequences that vary in their length distributions; this variance is a
    #' function of the sequencing library prepared, the starting DNA etc. This
    #' method is used to bin reads into uniform bins & assess the distribution
    #' of sequence lengths.
    #'
    #' @param qfilt - specifies how the quality information should be filtered
    #' - at the moment this only defines whether reporting is based on PASS
    #' or FAIL reads; would make more sense to have filtered by PASS / ALL?
    #' @param normalised - should the sequence collection be reported to
    #' normalise for the number of sequence bases sequenced or the number of
    #' sequence reads - TRUE by default to normalise for sequenced bases.
    #' @param cumulative defines whether cumulative sequence bases (reads) are
    #' reported per bin (FALSE by default).
    #' @param bins the number of sequence bins that should be prepared (20 by
    #' default)
    #' @param outliers defines the number of outliers (0.025 = 2.5%) that are
    #' excluded from the longest reads to prepare a richer distribution
    #' visulation - the plots can be bothered by the long tail of mini-whales.
    #'
    #' @return Angenieux 2D graph object
    read_length_bins = function(
      qfilt=TRUE, normalised=TRUE, cumulative=FALSE, bins=20, outliers=0.025) {
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
        rls <- rls %>%
          dplyr::group_by(.data[["barcode_arrangement"]]) %>%
          dplyr::arrange(desc(bin)) %>%
          dplyr::mutate(count=cumsum(count))
      }

      Angenieux$new("2D_count", rls)

    },


    #' @description
    #' bin the sequences in `seqsum` content into bins of quality
    #'
    #' The nanopore sequencing run is expected to return a collection of
    #' sequences that vary in their quality distributions; this variance is a
    #' function of the sequencing library prepared, the starting DNA etc. This
    #' method is used to bin reads into uniform quality bins to assess the
    #' overall quality of the run and to identify potential issues
    #'
    #' @param bins the number of sequence bins that should be prepared (20 by
    #' default)
    #' @param outliers defines the number of outliers (0 = 0%) that are
    #' excluded from the reads - should probably be deprecated for simplicity?
    #'
    #' @return Angenieux 2D graph object
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


    #' @description
    #' Prepare a SequencingSet object from a given barcode
    #'
    #' This method is used to subset the sequencing_summary information to
    #' focus on a single barcode for more detailed analysis.
    #'
    #' @param barcode the barcode that should be reported
    #'
    #' @return `SequencingSet` object containing information on barcode of
    #' interest.
    sequencingset = function(barcode) {
      if (!barcode %in% private$seqsum$barcode_arrangement) {
        stop(
          paste0(
            "barcode [",barcode,"] not present in sequencingsummary"))
      }
      return(SequencingSet$new(
        keycol = "passes_filtering",
        seqsum = private$seqsum %>%
          dplyr::filter(.data[["barcode_arrangement"]]==barcode)))
    },

    #' @description
    #' Prepare a TemporalSet object from a given barcode
    #'
    #' This method is used to subset the SequencingSummary information to
    #' focus on a single barcode for a more detailed analysis of content.
    #'
    #' @param barcode the barcode that should be reported
    #'
    #' @return `TemporalSet` object containing information on barcode of
    #' interest.
    temporalset = function(barcode) {
      if (!barcode %in% private$seqsum$barcode_arrangement) {
        stop(
          paste0(
            "barcode [",barcode,"] not present in sequencingsummary"))
      }
      return(TemporalSet$new(
        seqsum = private$seqsum %>%
          dplyr::filter(.data[["barcode_arrangement"]]==barcode)))
    }

  ),

  active = list (

    #' @field enumerate
    #' prepares a simple `2D Angenieux` enumeration of the provided dataset
    #' for quick visualisation of the dataset.
    enumerate = function(value) {
      if (missing(value)) {
        enumerated_counts <- private$seqsum %>%
          dplyr::group_by(
            passes_filtering, barcode_arrangement) %>%
          dplyr::summarize(count=dplyr::n(), .groups="drop")
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
