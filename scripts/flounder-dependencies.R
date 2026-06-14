flounder_cran_packages <- function() {
  c(
    "aws.s3",
    "cli",
    "dplyr",
    "emojifont",
    "ggplot2",
    "ggsci",
    "jsonlite",
    "knitr",
    "lubridate",
    "magick",
    "magrittr",
    "purrr",
    "R6",
    "RColorBrewer",
    "readr",
    "reshape2",
    "rlang",
    "rmarkdown",
    "stringdist",
    "stringr",
    "svglite",
    "testthat",
    "tibble",
    "tidyr",
    "tidyselect"
  )
}

flounder_bioc_packages <- function() {
  c(
    "BiocGenerics",
    "Biostrings",
    "GenomeInfoDb",
    "GenomicRanges",
    "IRanges",
    "Rsamtools",
    "ShortRead"
  )
}

flounder_parse_description_packages <- function(field) {
  if (is.na(field) || !nzchar(field)) {
    return(character())
  }

  packages <- trimws(unlist(strsplit(field, ",", fixed = TRUE)))
  packages <- sub("\\s*\\(.*\\)$", "", packages)
  packages[nzchar(packages)]
}

flounder_dependency_plan <- function(description_path = "DESCRIPTION") {
  description <- read.dcf(description_path)
  imports <- flounder_parse_description_packages(description[1, "Imports"])
  suggests <- flounder_parse_description_packages(description[1, "Suggests"])
  cran <- flounder_cran_packages()
  bioc <- flounder_bioc_packages()

  data.frame(
    package = unique(c(imports, suggests)),
    declared_in = ifelse(unique(c(imports, suggests)) %in% imports, "Imports", "Suggests"),
    support_channel = ifelse(unique(c(imports, suggests)) %in% bioc, "Bioconductor", "CRAN"),
    stringsAsFactors = FALSE
  )
}

flounder_expected_packages <- function() {
  c(flounder_cran_packages(), flounder_bioc_packages())
}

flounder_validate_dependency_plan <- function(description_path = "DESCRIPTION") {
  declared <- flounder_dependency_plan(description_path)$package
  missing_from_bootstrap <- setdiff(declared, flounder_expected_packages())

  if (length(missing_from_bootstrap) > 0) {
    stop(
      "Bootstrap list is missing DESCRIPTION dependencies: ",
      paste(missing_from_bootstrap, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
