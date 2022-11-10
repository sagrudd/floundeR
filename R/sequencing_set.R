
#' R6 Class for loading and analysing sequence sets
#'
#' @importFrom tidyr drop_na
#' @importFrom magrittr %>%
#' @importFrom dplyr across
#' @importFrom tidyselect all_of
#'
#' @export
SequencingSet <- R6::R6Class(
  inherit = FloundeR,
  classname = "SequencingSet",
  public = list(

    #' @description
    #' Initialise a new instance of the R6 Class `SequencingSet`
    #'
    #' @param keycol a pointer to the column of interest in the `seqsum` to
    #' direct parsing and exploration of the file.
    #' @param seqsum a tibble of sequencing summary information
    initialize = function(keycol, seqsum=NA) {
      if (is.character(keycol) & tibble::is_tibble(seqsum)) {
        if (! keycol %in% colnames(seqsum)) {
          stop("keycol must be a seqsum column name")
        }
        private$keycol <- keycol
        private$seqsum <- seqsum
      } else {

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
    #' method is used to bin reads into uniform bins to assess the distribution
    #' of sequence lengths.
    #'
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
      normalised=TRUE, cumulative=FALSE, bins=20, outliers=0.025) {
      bins <- self$bin_data(
        private$seqsum$sequence_length_template, bins=bins,
        outliers=outliers)
      private$seqsum["bin"] <- as.numeric(levels(bins))[bins]
      rls <- NULL
      if (normalised) {
        rls <- private$seqsum %>%
          dplyr::group_by(.dots=c(private$keycol, "bin")) %>%
          dplyr::summarize(count=sum(sequence_length_template), .groups="drop")
      } else {
        rls <- private$seqsum %>%
          dplyr::group_by(.dots=c(private$keycol, "bin")) %>%
          dplyr::summarize(count=dplyr::n(), .groups="drop")
      }
      # na may creep into the data ...
      rls <- tidyr::drop_na(rls)

      if (cumulative) {
        rls <- rls %>%
          dplyr::group_by(passes_filtering) %>%
          dplyr::arrange(desc(bin)) %>%
          dplyr::mutate(count=cumsum(count))
      }

      plot <- Angenieux$new("2D_count", rls)

      plot$add(
        list(
          AngenieuxDecoration$new("ggplot2", facet=xlab("Length of sequence (bases)")),
          AngenieuxDecoration$new("ggplot2",facet= scale_x_continuous(labels = self$nums_scale)),
          AngenieuxDecoration$new("vline", xintercept=self$N50, colour="gold", size=0.5),
          AngenieuxDecoration$new("vlegend", xintercept=self$N50, legend=" N50 ", size=3, colour="gold"),
          AngenieuxDecoration$new("vline", xintercept=self$mean, colour=plot$colourMap[3], size=0.5),
          AngenieuxDecoration$new("vlegend", xintercept=self$mean, legend=" mean ", colour=plot$colourMap[3], hjust=1, size=3),
          AngenieuxDecoration$new("ggplot2", facet=guides(shape = guide_legend(override.aes = list(size = 0.3)))),
          AngenieuxDecoration$new("ggplot2", facet=guides(color = guide_legend(override.aes = list(size = 0.3)))),
          AngenieuxDecoration$new("ggplot2", facet=theme(legend.title = element_text(size = 6), legend.text = element_text(size = 6)))
          )
      )


      if (normalised & cumulative) {
        plot$set_title("Plot showing normalised cumulative distribution of read lengths across sequence collection")
      } else if (normalised & !cumulative) {
        plot$set_title("Plot showing normalised distribution of read lengths across sequence collection")
      } else if (!normalised & cumulative) {
        plot$set_title("Plot showing cumulative distribution of sequenced reads by read length")
      } else if (!normalised & !cumulative) {
        plot$set_title("Plot showing distribution of sequenced reads by read length")
      }

      if (normalised) {
        plot$add(AngenieuxDecoration$new("ggplot2",facet = scale_y_continuous(labels = self$nums_scale)))
        plot$add(AngenieuxDecoration$new("ggplot2",facet = ylab("Sequence bases")))
      } else {
        plot$add(AngenieuxDecoration$new("ggplot2",facet = ylab("Sequence reads")))
      }

      return(plot)
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
    #' excluded from the reads - should probably be deprecated for simplicity??
    #'
    #' @return Angenieux 2D graph object
    quality_bins = function(bins=20, outliers=0) {
      bins <- self$bin_data(
        private$seqsum$mean_qscore_template, bins=bins,
        outliers=outliers)
      private$seqsum["bin"] <- as.numeric(levels(bins))[bins]
      rls <- private$seqsum %>%
        dplyr::group_by(.dots=c(private$keycol, "bin")) %>%
        dplyr::summarize(count=dplyr::n(), .groups="drop")
      # na may creep into the data ...
      rls <- tidyr::drop_na(rls)
      rls$bin <-  factor(rls$bin, sort(unique(rls$bin)))
      Angenieux$new("2D_count", rls)
    }
  ),

  active = list (

    #' @field enumerate
    #' prepares a simple `1D Angenieux` enumeration of the provided dataset
    #' for quick visualisation of the dataset.
    enumerate = function(value) {
      if (missing(value)) {
        enumerated_counts <- private$seqsum %>%
          dplyr::group_by(dplyr::across(tidyselect::all_of(private$keycol))) %>%
          dplyr::summarize(count=dplyr::n())
        plot <- Angenieux$new("1D_count", enumerated_counts)
        plot$set_title(
          stringr::str_interp(
            "Plot showing count information for reads that ${private$keycol}"))
        plot$add(
          list(
            AngenieuxDecoration$new("ggplot2",facet = ylab("Sequence reads")),
            AngenieuxDecoration$new("ggplot2",facet= scale_y_continuous(labels = self$nums_scale))
            )
        )

        return(plot)
      }
    },

    #' @field N50
    #' Calculate and return the N50 value for passed quality sequence reads in
    #' the current `SequencingSet` object
    N50 = function(value) {
      if (missing(value)) {
        len_sorted <- private$seqsum %>%
          dplyr::filter(passes_filtering == TRUE) %>%
          dplyr::select(sequence_length_template) %>%
          dplyr::arrange(desc(sequence_length_template))
        len_sorted[cumsum(len_sorted) >= sum(len_sorted)*0.5][1]
      }
    },

    #' @field mean
    #' Calculate and return the mean sequence length for passed quality reads in
    #' the `SequencingSet` object
    mean = function(value) {
      if (missing(value)) {
        as.numeric(
          private$seqsum %>%
            dplyr::filter(passes_filtering == TRUE) %>%
            dplyr::select(sequence_length_template) %>%
            dplyr::summarize(mean(sequence_length_template)))
      }
    }

  ),

  private = list(
    keycol = NULL,
    seqsum = NULL
  ),
)





#' Prepare a `SequencingSet` object from provided class
#'
#' This method is intended to streamline the usage of floundeR in more coherent
#' workflows through the enabling of more magrittr-like piped commands. This
#' method will consume floundeR specified classes and return a SequencingSet
#' object where possible.
#'
#' @param r6_object is the floundeR R6 object to extract SequencingSet from
#'
#' @return Flounder SequencingSet object
#'
#' @importFrom magrittr %>%
#' @export %>%
#'
#' @examples
#' SequencingSummary$new(flnDr("sequencing_summary.txt.bz2")) %>%
#'     to_sequencing_set()
#'
#' @export
to_sequencing_set = function(r6_object) {
  if (class(r6_object)[1] == "SequencingSummary") {
    return(r6_object$sequencingset)
  } else if (class(r6_object)[1] == "Fastq") {
    return(r6_object$sequencingset)
  }  else if (class(r6_object)[1] == "Fasta") {
    return(r6_object$sequencingset)
  }  else {
    stop("Unable to prepare flowcell object from provided input")
  }

}
