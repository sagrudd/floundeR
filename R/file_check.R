
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
