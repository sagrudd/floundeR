

#' Sanity check a file and ensure that it is a valid FAST5 file
#'
#' Simple method to check whether the provided file path corresponds to a valid
#' FAST5 file (or not). This method operates only on atomic files.
#'
#' @param path The file.path to check
#' @return boolean
#'
#' @examples
#' singleFast5 <- system.file("extdata", "single.fast5", package="floundeR")
#' file.isFast5(singleFast5)
#'
#' @export
file.isFast5 <- function(path) {
    f5content = file.isHdf5(path)
    if (is.logical(f5content)) {
        return(f5content)
    }
    if (file.isSingleFast5(f5content=f5content, silent=TRUE)) {
        return(TRUE)
    } else if (file.isMultiFast5(f5content=f5content, silent=TRUE)) {
        return(TRUE)
    }
    message(paste0("path [",path,"] provides HDF but not nanopore FAST5"))
    return(FALSE)
}

#' @importFrom utils file_test
#' @importFrom rhdf5 h5ls
file.isHdf5 <- function(path) {
    if (!file.exists(path)) {
        message(paste0("path [",path,"] does not exist"))
        return(FALSE)
    } else if (!utils::file_test("-f", path)) {
        message(paste0("path [",path,"] is a directory"))
        return(FALSE)
    }
    f5content = tryCatch(
        rhdf5::h5ls(path),
        error=function(err){return(NULL)},
        finally={}
    )
    if (is.null(f5content)) {
        message(paste0("file [",path,"] cannot be parsed with rhdf5"))
        return(FALSE)
    } else {
        return(f5content)
    }

}


getH5list <- function(path=NULL, f5content=NULL) {
    if (!is.null(path) & is.null(f5content)) {
        f5content = file.isHdf5(path)
        if (is.logical(f5content)) {
            return(NULL)
        }
    }
    return(f5content)
}

#' Check is a provided file corresponds to a single entry FAST5 file
#'
#' Simple method to check whether the provided file path corresponds to a valid
#' FAST5 file and if this file corresponds to a single FAST5 format.
#'
#' @param path The file.path to check
#' @param f5content (optional) preparsed index from h5parser
#' @param silent should messages be hidden? (FALSE)
#' @return boolean
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#'
#' @examples
#' singleFast5 <- system.file("extdata", "single.fast5", package="floundeR")
#' file.isSingleFast5(singleFast5)
#' multiFast5 <- system.file("extdata", "multi.fast5", package="floundeR")
#' file.isSingleFast5(multiFast5)
#'
#' @export
file.isSingleFast5 <- function(path=NULL, f5content=NULL, silent=FALSE) {
    f5content = getH5list(path, f5content)
    if (is.null(f5content)) {
        message("something has gone awry")
        return(FALSE)
    }

    if (all(
        c("Analyses", "Raw", "UniqueGlobalKey") %in%
        dplyr::filter(f5content, .data$group=="/")[["name"]])) {
        return(TRUE)
    }
    if (!silent) message("This appears to be HDF5 but not single FAST5")
    return(FALSE)
}

#' Check is a provided file corresponds to a multiple entry FAST5 file
#'
#' Simple method to check whether the provided file path corresponds to a valid
#' FAST5 file and if this file corresponds to a multi FAST5 format.
#'
#' @param path The file.path to check
#' @param f5content (optional) preparsed index from h5parser
#' @param silent should messages be hidden? (FALSE)
#' @return boolean
#'
#' @seealso \code{\link{file.isSingleFast5} \link{file.isFast5}}
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @examples
#' singleFast5 <- system.file("extdata", "single.fast5", package="floundeR")
#' file.isMultiFast5(singleFast5)
#' multiFast5 <- system.file("extdata", "multi.fast5", package="floundeR")
#' file.isMultiFast5(multiFast5)
#'
#' @export
file.isMultiFast5 <- function(path=NULL, f5content=NULL, silent=FALSE) {
    f5content = getH5list(path, f5content)
    if (is.null(f5content)) {
        message("something has gone awry")
        return(FALSE)
    }
    itemA = dplyr::filter(f5content, .data$group=="/")[["name"]][1]
    if (all(
        c("Raw", "channel_id", "context_tags", "tracking_id") %in%
        dplyr::filter(f5content, .data$group==paste0("/",itemA))[["name"]])) {
        return(TRUE)
    }
    if (!silent) message("This appears to be HDF5 but not multi FAST5")
    return(FALSE)

}


#' Extract the sequencing platform identified within a FAST5 file
#'
#' This method parses a FAST5 sequence file and extracts the information on
#' specified sequencing platform - this can be used to identify e.g. gridion or
#' promethion based information.
#'
#' @param path The file.path to check
#' @param atomic whether a single entity (multi files only) should be considered
#'     (TRUE)
#' @return vector of platform(s)
#'
#' @examples
#' promFast5 <- system.file("extdata", "prom.fast5", package="floundeR")
#' fast5.get_platform(promFast5)
#'
#' @export
fast5.get_platform <- function(path, atomic=TRUE) {
    anno <- fast5.get_info(path, atomic)
    return(anno[,"device_type"])
}


#' Extract the flowcell_id identified within a FAST5 file
#'
#' This method parses a FAST5 sequence file and extracts the information on
#' specified flowcell_id - this can be used to track flowcells used within a
#' laboratory
#'
#' @param path The file.path to check
#' @param atomic whether a single entity (multi files only) should be considered
#'     (TRUE)
#' @return vector of flowcell_id(s)
#'
#' @examples
#' multiFast5 <- system.file("extdata", "multi.fast5", package="floundeR")
#' fast5.get_flowcell_id(multiFast5)
#'
#' @export
fast5.get_flowcell_id <- function(path, atomic=TRUE) {
    anno <- fast5.get_info(path, atomic)
    return(anno[,"flow_cell_id"])
}


#' Extract the experiment start time identified within a FAST5 file
#'
#' This method parses a FAST5 sequence file and extracts the information on
#' specified experiment start time - this can be used in palaeogenomics of such
#' ancient datasets
#'
#' @param path The file.path to check
#' @param atomic whether a single entity (multi files only) should be considered
#'     (TRUE)
#' @return vector of experiment_start_time(s)
#'
#' @examples
#' multiFast5 <- system.file("extdata", "multi.fast5", package="floundeR")
#' fast5.get_exp_start_time(multiFast5)
#'
#' @export
fast5.get_exp_start_time <- function(path, atomic=TRUE) {
    anno <- fast5.get_info(path, atomic)
    return(lubridate::ymd_hms(anno$exp_start_time))
}


#' Extract the number of reads identified within a FAST5 file
#'
#' This method parses a FAST5 sequence file and extracts the information on
#' the read count - this is a simple count only
#'
#' @param path The file.path to check
#' @return integer count of reads
#'
#' @examples
#' multiFast5 <- system.file("extdata", "multi.fast5", package="floundeR")
#' fast5.get_read_count(multiFast5)
#'
#' @export
fast5.get_read_count <- function(path) {
    f5content = getH5list(path)
    if (is.null(f5content)) {
        message("something has gone awry")
        return(0)
    }
    if (file.isSingleFast5(f5content=f5content, silent=TRUE)) {
        return(1)
    } else if (file.isMultiFast5(f5content=f5content, silent=TRUE)) {
        named_sequences = dplyr::filter(f5content, .data$group=="/")[["name"]]
        return(length(named_sequences))
    }
    return(0)
}



#' Extract the flowcell type identified within a FAST5 file
#'
#' This method parses a FAST5 sequence file and extracts the information on
#' the flowcell type - this is based on user entry at start of run.
#'
#' @param path The file.path to check
#' @param atomic whether a single entity (multi files only) should be considered
#'     (TRUE)
#' @return vector of flowcell_types(s)
#'
#' @examples
#' multiFast5 <- system.file("extdata", "multi.fast5", package="floundeR")
#' fast5.get_flowcell_type(multiFast5)
#'
#' @export
fast5.get_flowcell_type <- function(path, atomic=TRUE) {
    anno <- fast5.get_info(path, atomic)
    return(anno[,"flowcell_type"])
}


#' Extract the sequencing kit identified within a FAST5 file
#'
#' This method parses a FAST5 sequence file and extracts the information on
#' the sequencing_kit used - this is based on user entry at start of run.
#'
#' @param path The file.path to check
#' @param atomic whether a single entity (multi files only) should be considered
#'     (TRUE)
#' @return vector of sequencing_kit(s)
#'
#' @examples
#' multiFast5 <- system.file("extdata", "multi.fast5", package="floundeR")
#' fast5.get_sequencing_kit(multiFast5)
#'
#' @export
fast5.get_sequencing_kit <- function(path, atomic=TRUE) {
    anno <- fast5.get_info(path, atomic)
    return(anno[,"sequencing_kit"])
}



#' @importFrom rhdf5 h5closeAll
#' @importFrom rhdf5 h5readAttributes
fast5.get_attributes <- function(key, path) {
    tracking_attrs <- rhdf5::h5readAttributes(
        path, paste0("/", key, "/tracking_id"))
    rhdf5::h5closeAll()
    context_attrs <- rhdf5::h5readAttributes(
        path, paste0("/", key, "/context_tags"))
    rhdf5::h5closeAll()
    return(data.frame(
        device_type=tracking_attrs$device_type,
        flow_cell_id=tracking_attrs$flow_cell_id,
        exp_start_time=tracking_attrs$exp_start_time,
        experiment_type=context_attrs$experiment_type,
        flowcell_type=context_attrs$flowcell_type,
        sequencing_kit=context_attrs$sequencing_kit))
}


#' Extract salient experimental information from a FAST5 file
#'
#' The FAST5 sequence file contains information relating to the sequencing
#' platform, flowcell and library preparation kits used in a study. This method
#' pulls these information into a data.frame or factor.
#'
#' @param path The file.path to check
#' @param atomic whether a single entity (multi files only) should be considered
#'     (FALSE)
#' @return data.frame
#'
#' @seealso \code{\link{file.isSingleFast5} \link{file.isFast5}}
#' @importFrom purrr map
#' @importFrom dplyr filter
#' @importFrom stats runif
#'
#' @examples
#' promFast5 <- system.file("extdata", "prom.fast5", package="floundeR")
#' fast5.get_info(promFast5)
#'
#' @export
fast5.get_info <- function(path, atomic=FALSE) {
    f5content = getH5list(path)
    if (is.null(f5content)) {
        message("something has gone awry")
        return(NULL)
    }
    lookup_key = NULL
    if (file.isSingleFast5(f5content=f5content, silent=TRUE)) {
        lookup_key = "UniqueGlobalKey"
    } else if (file.isMultiFast5(f5content=f5content, silent=TRUE)) {
        lookup_key = dplyr::filter(f5content, .data$group=="/")[["name"]]
        if (atomic) {
            pick = floor(runif(1, min=1, max=length(lookup_key)+1))
            lookup_key = lookup_key[pick]
        }
    }
    if (is.null(lookup_key)) {
        message(paste0("path [",path,"] provides HDF but not nanopore FAST5"))
        return(NULL)
    }
    l <- purrr::map_df(lookup_key, ~fast5.get_attributes(.x, path=path))
    return(l)

}


