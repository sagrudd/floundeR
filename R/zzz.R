
.onAttach <- function(libname, pkgname) {
    packageStartupMessage(paste0("floundeR v",
                                 utils::packageVersion("floundeR")))
    invisible()
}


.onLoad <- function(libname, pkgname){
    library(magrittr)
    invisible()
}
