
#' @importFrom utils file_test
.check_path = function(key, file.path) {
    if (is.na(file.path)) {
        stop(paste0(key, " file must be defined"))
    } else if (!is.character(file.path) ||
               length(file.path) != 1) {
        stop(paste0(key," requires a single [file.path] as input"))
    } else if (!file.exists(file.path)) {
        stop(paste0("path [",file.path,"] does not exist"))
    } else if (!utils::file_test("-f", file.path)) {
        stop(paste0("path [",file.path,
                    "] is a directory - file reqd"))
    }
}


#' extract a floundeR packaged file
#'
#' A collection of files are distributed with the floundeR package; this
#' accessory method is aimed to abbreviate the technical code for the
#' identification of these files and to clean the vignettes and code examples.
#'
#' @param file a known filename (or its prefix)
#'
#' @return a file
#'
#' @export
flnDr = function(file) {
    datadir = file.path(system.file(package="floundeR"), "extdata")
    return(
        file.path(datadir, list.files(datadir, pattern=file)[1])
    )
}
