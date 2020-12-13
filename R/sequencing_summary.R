
#' R6 Class for loading and analysing nanopore sequencing_summary files
#'
#' @description
#' This class aims to simplify the handling and exploration of
#' sequencing_summary files and provides simple methods for accessing
#' information that can be used to assess the performance of a run.
#'
#' @import R6
#' @importFrom reshape2 acast
#'
#' @export
SequencingSummary <- R6::R6Class(
    "SequencingSummary",
    public = list(
        #' @field sequencing_summary_file the file.path
        #' to the query FAST5 file
        sequencing_summary_file = NULL,
        #' @field barcoding_summary_file the file.path
        #' to the query FAST5 file
        barcoding_summary_file = NULL,

        #' @description
        #' Creates a new SequencingSummary object. This
        #' initialisation method performs other sanity checking
        #' of the defined file(s) to ensure that it is indeed
        #' parseable and creates the required data structures.
        #'
        #' @param sequencing_summary_file The source
        #' sequencing_summary file.
        #' @param barcoding_summary_file The source
        #' barcoding_summary_file file.
        #' @return A new `Fast5` object.
        #'
        #' @examples
        #' sequencing_summary <- flnDr("sequencing_summary.txt.bz2")
        #' barcodes_summary <- flnDr("barcoding_summary.txt.bz2")
        #' seqsum <- SequencingSummary$new(sequencing_summary, barcodes_summary)
        initialize = function(
            sequencing_summary_file,
            barcoding_summary_file = NA) {

            # first pass QC - let's ensure that this is a single FILE
            .check_path(
                "sequencing_summary_file", sequencing_summary_file)
            self$sequencing_summary_file = sequencing_summary_file

            if (!is.na(barcoding_summary_file)) {
                .check_path(
                    "barcoding_summary_file", barcoding_summary_file)
                self$barcoding_summary_file = barcoding_summary_file
            }
            # and import the TSV data (cleanly)
            if (!private$.parse_seqsum()) {
                stop("Failed to import the summary files provided ...")
            }
        },

        #' @description
        #' Have a guess at the most likely flowcell platform used
        #'
        #' The sequencing summary file contains no information on the sequencing
        #' device or flowcell used. For the preparation of channel density maps
        #' it is worth considering which flowcell type is most likely to have
        #' been used - this can be guessed on the number of channels described
        #' within the data
        #'
        #' @return character
        #'
        #' @examples
        #' seqsum <- SequencingSummary$new(flnDr("sequencing_summary.txt.bz2"))
        #' seqsum.get_flowcell_platform()
        get_flowcell_platform = function() {
            if (private$seqsum_channel_max < 130) {
                return("Flongle")
            } else if (private$seqsum_channel_max > 1000) {
                return("PromethION")
            }
            return("MinION")

        },


        #' @description
        #' produce channelMap for the predicted flowcell type for spatial plots
        #'
        #' prepares a matrix of X, Y coordinates and
        #' the corresponding readcount information for the type of flowcell
        #' predicted by `get_flowcell_platform`
        #'
        #' @return matrix with counts sorted by spatial rows and columns
        #'
        #' @examples
        #' seqsum <- SequencingSummary$new(flnDr("sequencing_summary.txt.bz2"))
        #' seqsum.get_fc_density_data()
        get_fc_density_data = function() {
            channelMap <- private$.get_channel_counts()
            channelMapMatrix <-
                reshape2::acast(channelMap, col ~ row, value.var = "count")
            return(XYDensity$new(channelMapMatrix))
        }


    ),

    private = list(

        seqsum = NULL,
        seqsum_channel_max = NULL,
        # read_id / c has been removed - this is a big character - ?value
        select_columns = c(
            "channel", "start_time", "duration",
            "passes_filtering", "sequence_length_template",
            "mean_qscore_template", "barcode_arrangement"),

        select_column_types = c(
            "i", "d", "d", "l", "d", "d", "c"
        ),

        .parse_seqsum = function() {

            mini_table = readr::read_tsv(
                file=self$sequencing_summary_file, n_max=10, col_types = readr::cols())
            import_cols <- names(mini_table)[
                na.omit(match(private$select_columns, names(mini_table)))]
            ctypes = private$select_column_types[
                match(import_cols, private$select_columns)]
            colt = paste0(
                "readr::cols_only(",
                paste(paste(
                    import_cols, "=\"", ctypes, "\"",sep=""), collapse=", "),
                ")")

            private$seqsum <- readr::read_tsv(
                file=self$sequencing_summary_file,
                col_types=eval(parse(text=colt)))

            # set some key values
            private$seqsum_channel_max = max(private$seqsum$channel)
            return(TRUE)
        },

        .get_minion_channel_map = function() {
            # build the map for R9.4.1 flowcell, as a long-form dataframe
            blockCalc <- function(i) {
                m <- matrix(seq(i, i + 63, by = 1), ncol = 8, byrow = TRUE)
                cbind(m[seq(5, 8, by = 1), ], m[seq(4), rev(seq(8))])
            }
            layout <- do.call(rbind, lapply(
                c(1, 449, 385, 321, 257, 193, 129, 65), blockCalc))
            # transpose the layout for cleaner presentation ...
            layout <- t(layout)
            channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)), which(
                layout == as.vector(layout), arr.ind = TRUE)))
            return(channelMap)
        },

        .get_flongle_channel_map = function() {
            layout <- matrix(c(seq(1, 12), 0, seq(13, 24), 0, seq(25, 114), 0,
                               seq(115, 126), 0), ncol = 13, byrow = TRUE)
            layout <- layout[rev(seq(10)), ]
            channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)),
                                              which(layout == as.vector(layout), arr.ind = TRUE)))
            return(channelMap)
        },

        .get_promethion_channel_map = function() {
            chunk <- function(i) {
                m <- matrix(seq_len(250), ncol=10, byrow=TRUE)
                m + i
            }
            layout <- do.call(cbind, lapply(seq(from=0, to=2750, by=250), chunk))
            channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)),
                                              which(layout == as.vector(layout), arr.ind = TRUE)))
            return(channelMap)
        },

        .get_channel_map = function() {
            platform = self$get_flowcell_platform()
            if (platform == "MinION") {
                return(private$.get_minion_channel_map())
            } else if (platform == "Flongle") {
                return(private$.get_flongle_channel_map())
            } else if (platform == "PromethION") {
                return(private$.get_promethion_channel_map())
            }
        },

        .get_channel_counts = function() {
            channelMap <- private$.get_channel_map()
            channelCounts <-
                as.data.frame(matrix(rep(0, max(channelMap$channel)), ncol = 1))
            channelCountRaw <-
                as.data.frame(
                    table(unlist(private$seqsum[, "channel"])), row.names = 1)
            channelCounts[row.names(channelCountRaw), ] <- channelCountRaw[, 1]
            channelMap <- merge(channelMap, channelCounts, by.x = "channel", by.y = 0)
            colnames(channelMap)[4] <- "count"
            return(channelMap)
        }

    )
)


XYDensity <- R6::R6Class(
    "XYDensity",
    public = list(
        matrix = NULL,

        initialize = function(matrix) {
            if (!is.matrix(matrix)) {
                stop("this requires a matrix")
            }
            self$matrix <- matrix
        },

        print = function(...) {
            cat(paste0("floundeR::XYDensity[",
                       ncol(self$matrix),", ",nrow(self$matrix),"]"))
        },

        data = function() {
            return(self$matrix)
        },

        plot = function() {
            print("a plot")
        }

    )
)


