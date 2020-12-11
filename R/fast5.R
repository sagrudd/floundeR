#' R6 Class for loading and analysing nanopore FAST5 sequence files
#'
#' @description
#' This class aims to simplify the handling and exploration of FAST5 files and
#' provides simple methods for accessing information that can be used to
#' identify the sequencing platform, flowcell used, sequencing kit and other
#' metainformation that includes whether they are single- or- multi format
#' files.
#'
#' @import R6
#' @importFrom utils file_test
#' @importFrom rhdf5 h5ls
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom rhdf5 h5closeAll
#' @importFrom rhdf5 h5readAttributes
#' @importFrom rhdf5 h5closeAll
#' @importFrom stats runif
#' @importFrom purrr map
#' @importFrom lubridate ymd_hms
#'
#' @export
Fast5 <- R6::R6Class(
    "Fast5",
    public = list(
        #' @field fast5_file the file.path to the query FAST5 file
        fast5_file = NULL,

        #' @description
        #' Creates a new Fast5 object. This initialisation method performs other
        #' sanity checking of the defined file to ensure that it is indeed a
        #' file, a single file, is in HDF5 format and is either single- or
        #' multi-FAST5 compliant
        #'
        #' @param fast5_file The source FAST5 file.
        #' @return A new `Fast5` object.
        initialize = function(fast5_file = NA) {
            # first pass QC - let's ensure that this is a single FILE
            if (!is.character(fast5_file) || length(fast5_file) != 1) {
                stop(paste0("FAST5 requires a single [file.path] as input"))
            } else if (!file.exists(fast5_file)) {
                stop(paste0("path [",fast5_file,"] does not exist"))
            } else if (!utils::file_test("-f", fast5_file)) {
                stop(paste0("path [",fast5_file,"] is a directory - file reqd"))
            }

            # second pass QC- ensure file is actually a parseable FAST5 file
            private$f5content = tryCatch(
                rhdf5::h5ls(fast5_file),
                error=function(err){return(NULL)},
                finally={}
            )
            if (is.null(private$f5content)) {
                stop(paste0("file [",fast5_file,
                            "] cannot be parsed with rhdf5"))
            }

            # is this a single or multi FAST5 file?

            if (private$.is_singleF5()) {
                self$fast5_file <- fast5_file
                private$singleF5 = TRUE
                private$atomic_anno = private$.get_info(atomic=TRUE)
            } else if (private$.is_multiF5()) {
                self$fast5_file <- fast5_file
                private$multiF5 = TRUE
                private$atomic_anno = private$.get_info(atomic=TRUE)
            } else {
                stop(paste0("[",fast5_file,
                            "] appears to be HDF5 but not FAST5"))
            }
        },

        #' @description
        #' Check if a provided file corresponds to a single entry FAST5 file
        #'
        #' Simple method to check whether the provided file path corresponds to
        #' a valid FAST5 file and if this file corresponds to a single FAST5
        #' format.
        #'
        #' @return boolean
        #'
        #' @examples
        #' promFast5 <- system.file("extdata", "prom.fast5", package="floundeR")
        #' promF5 <- Fast5$new(fast5_file=promFast5)
        #' promF5$is_single_fast5()
        is_single_fast5 = function() {
            return(private$singleF5)
        },

        #' @description
        #' Check if a provided file corresponds to a multiple entry FAST5 file
        #'
        #' Simple method to check whether the provided file path corresponds to
        #' a valid FAST5 file and if this file corresponds to a multi FAST5
        #' format.
        #'
        #' @return boolean
        #'
        #' @examples
        #' promFast5 <- system.file("extdata", "prom.fast5", package="floundeR")
        #' promF5 <- Fast5$new(fast5_file=promFast5)
        #' promF5$is_multi_fast5()
        is_multi_fast5 = function() {
            return(private$multiF5)
        },

        #' @description
        #' Extract the sequencing platform identified within a FAST5 file
        #'
        #' This method parses a FAST5 sequence file and extracts the information
        #' on specified sequencing platform - this can be used to identify e.g.
        #' gridion or promethion based information.
        #'
        #' @return vector of platform
        #'
        #' @examples
        #' promFast5 <- system.file("extdata", "prom.fast5", package="floundeR")
        #' promF5 <- Fast5$new(fast5_file=promFast5)
        #' promF5$get_platform()
        get_platform = function() {
            return(private$atomic_anno[,"device_type"])
        },

        #' @description
        #' Extract the flowcell_id identified within a FAST5 file
        #'
        #' This method parses a FAST5 sequence file and extracts the information
        #' on specified flowcell_id - this can be used to track flowcells used
        #' within a laboratory
        #'
        #' @return vector of flowcell_id
        #'
        #' @examples
        #' promFast5 <- system.file("extdata", "prom.fast5", package="floundeR")
        #' promF5 <- Fast5$new(fast5_file=promFast5)
        #' promF5$get_flowcell_id()
        get_flowcell_id = function() {
            return(private$atomic_anno[,"flow_cell_id"])
        },

        #' @description
        #' Extract the experiment start time identified within a FAST5 file
        #'
        #' This method parses a FAST5 sequence file and extracts the information
        #' on specified experiment start time - this can be used in
        #' palaeogenomics of such ancient datasets. Please note that the sampled
        #' start time is transformed using `lubridate` for presentation.
        #'
        #' @return vector of experiment_start_time
        #'
        #' @examples
        #' promFast5 <- system.file("extdata", "prom.fast5", package="floundeR")
        #' promF5 <- Fast5$new(fast5_file=promFast5)
        #' promF5$get_exp_start_time()
        get_exp_start_time = function() {
            return(lubridate::ymd_hms(private$atomic_anno[,"exp_start_time"]))
        },


        #' @description
        #' Extract the number of reads identified within a FAST5 file
        #'
        #' This method parses a FAST5 sequence file and extracts the information
        #' on the read count - this is a simple count only
        #'
        #' @return integer count of reads
        #'
        #' @examples
        #' promFast5 <- system.file("extdata", "prom.fast5", package="floundeR")
        #' promF5 <- Fast5$new(fast5_file=promFast5)
        #' promF5$get_read_count()
        get_read_count = function() {
            if (private$.is_singleF5()) {
                return(1)
            } else if (private$.is_multiF5()) {
                return(length(private$lookup_key))
            }
        },


        #' @description
        #' Extract the flowcell type identified within a FAST5 file
        #'
        #' This method parses a FAST5 sequence file and extracts the information
        #' on the flowcell type - this is based on user entry at start of run.
        #'
        #' @return vector of flowcell_types
        #'
        #' @examples
        #' promFast5 <- system.file("extdata", "prom.fast5", package="floundeR")
        #' promF5 <- Fast5$new(fast5_file=promFast5)
        #' promF5$get_flowcell_type()
        get_flowcell_type = function() {
            return(private$atomic_anno[,"flowcell_type"])
        },


        #' @description
        #' Extract the sequencing kit identified within a FAST5 file
        #'
        #' This method parses a FAST5 sequence file and extracts the information
        #' on the sequencing_kit used - this is based on user entry at start of
        #' run.
        #'
        #' @return vector of sequencing_kit
        #'
        #' @examples
        #' promFast5 <- system.file("extdata", "prom.fast5", package="floundeR")
        #' promF5 <- Fast5$new(fast5_file=promFast5)
        #' promF5$get_sequencing_kit()
        get_sequencing_kit = function() {
            return(private$atomic_anno[,"sequencing_kit"])
        },

        #' @description
        #' Extract salient experimental information from a FAST5 file
        #'
        #' The FAST5 sequence file contains information relating to the
        #' sequencing platform, flowcell and library preparation kits used in a
        #' study. This method pulls these information into a data.frame.
        #'
        #' @param atomic whether a single entity (multi files only) should be
        #' considered (FALSE)
        #' @return data.frame
        #'
        #' @examples
        #' promFast5 <- system.file("extdata", "prom.fast5", package="floundeR")
        #' promF5 <- Fast5$new(fast5_file=promFast5)
        #' promF5$get_info(atomic=FALSE)
        get_info = function(atomic=TRUE) {
            if (atomic) {
                return(private$atomic_anno)
            } else {
                return(private$.get_info(atomic))
            }
        },

        #' @description
        #' print method to override the standard R6 printing
        #'
        #' @param ... stuff passed onwards
        print = function(...) {
            cat("floundeR::Fast5<R6>Object\n")
            print(as_tibble(self$get_info()))
        }

    ),

    private = list(
        f5content = NULL,
        singleF5 = FALSE,
        multiF5 = FALSE,
        atomic_anno = NULL,
        pick = NULL,
        lookup_key = NULL,


        .is_singleF5 = function() {
            return(all(
                c("Analyses", "Raw", "UniqueGlobalKey") %in%
                    dplyr::filter(private$f5content,
                                  .data$group=="/")[["name"]]))
        },


        .is_multiF5 = function() {
            private$lookup_key = dplyr::filter(
                private$f5content, .data$group=="/")[["name"]]
            private$pick = private$lookup_key[
                floor(runif(1, min=1, max=length(private$lookup_key)+1))]

            return(all(
                c("Raw", "channel_id", "context_tags", "tracking_id") %in%
                    dplyr::filter(
                        private$f5content,
                        .data$group==paste0("/",private$pick))[["name"]]))
        },


        .get_attributes = function(key) {
            tracking_attrs <- rhdf5::h5readAttributes(
                self$fast5_file, paste0("/", key, "/tracking_id"))
            rhdf5::h5closeAll()
            context_attrs <- rhdf5::h5readAttributes(
                self$fast5_file, paste0("/", key, "/context_tags"))
            rhdf5::h5closeAll()
            return(data.frame(
                device_type=tracking_attrs$device_type,
                flow_cell_id=tracking_attrs$flow_cell_id,
                exp_start_time=tracking_attrs$exp_start_time,
                experiment_type=context_attrs$experiment_type,
                flowcell_type=context_attrs$flowcell_type,
                sequencing_kit=context_attrs$sequencing_kit))
        },

        .get_info = function(atomic=FALSE) {
            if (private$singleF5) {
                return(private$.get_attributes("UniqueGlobalKey"))
            } else if (private$multiF5 && atomic) {
                return(private$.get_attributes(private$pick))
            } else if (private$multiF5 && !atomic) {
                return(
                    purrr::map_df(
                        private$lookup_key,
                        ~private$.get_attributes(.x)))
            }
        }

    )

)




