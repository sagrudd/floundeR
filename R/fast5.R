

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
    if (file.isSingleFast5(f5content=f5content)) {
        return(TRUE)
    } else if (file.isMultiFast5(f5content=f5content)) {
        return(TRUE)
    }
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
            return(f5content)
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
file.isSingleFast5 <- function(path=NULL, f5content=NULL) {
    f5content = getH5list(path, f5content)
    if (is.null(f5content)) {
        message("something has gone awry")
        return(FALSE)
    }

    if (all(
        c("Analyses", "Raw", "UniqueGlobalKey") %in%
        dplyr::filter(f5content, .data$group=="/")[["name"]])) {
        message("This appears to be a single FAST5 file")
        return(TRUE)
    }
    message("This appears to be HDF5 but not single FAST5")
    return(FALSE)
}

#' Check is a provided file corresponds to a multiple entry FAST5 file
#'
#' Simple method to check whether the provided file path corresponds to a valid
#' FAST5 file and if this file corresponds to a multi FAST5 format.
#'
#' @param path The file.path to check
#' @param f5content (optional) preparsed index from h5parser
#' @return boolean
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
file.isMultiFast5 <- function(path=NULL, f5content=NULL) {
    f5content = getH5list(path, f5content)
    if (is.null(f5content)) {
        message("something has gone awry")
        return(FALSE)
    }
    itemA = dplyr::filter(f5content, .data$group=="/")[["name"]][1]
    if (all(
        c("Analyses", "Raw", "channel_id", "context_tags", "tracking_id") %in%
        dplyr::filter(f5content, .data$group==paste0("/",itemA))[["name"]])) {
        message("This appears to be a multi FAST5 file")
        return(TRUE)
    }
    message("This appears to be HDF5 but not multi FAST5")
    return(FALSE)

}


fast5info <- function(path) {

}
