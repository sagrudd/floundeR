# TB_reference = flnDr("NC_000962")
# GenbankGenome$new(TB_reference)

#' @export
GenbankGenome <- R6::R6Class(
    inherit = FloundeR,
    classname = "GenbankGenome",
    public = list(

        initialize = function(gb_file) {
            private$process_file(gb_file)
        },



        get_cds = function() {
            return(
                GenomicRanges::makeGRangesFromDataFrame(
                    private$gb_cds, keep.extra.columns = TRUE))
        },


        get_sequence = function() {

        }



    ),

    active = list (

        accession = function(set=NULL) {
            if (!is.null(set)) {
                cli::cli_alert(stringr::str_interp("setting accession [${set}]"))
                private$gb_accession <- set
            }
            return(private$gb_accession)
        },

        version = function(set=NULL) {
            if (!is.null(set)) {
                cli::cli_alert(stringr::str_interp("setting version [${set}]"))
                private$gb_version <- set
            }
            return(private$gb_version)
        },

        definition = function(set=NULL) {
            if (!is.null(set)) {
                cli::cli_alert(stringr::str_interp("setting definition [${set}]"))
                private$gb_definition <- set
            }
            return(private$gb_definition)
        }

    ),

    private = list(
        conn = NA,
        current_tag = NA,
        tag_string = NULL,
        key_tags = c("DEFINITION", "ACCESSION", "VERSION", "FEATURES", "ORIGIN"),
        debug_counter = 0,
        gb_cds = data.frame(
            gene=character(0),
            chr=character(0),
            strand=character(0),
            start=integer(0),
            end=integer(0),
            translation=character(0)),
        gb_accession = "undefined",
        gb_version = "undefined",
        gb_definition = "undefined",

        process_file = function(gb_file) {
            line_count <- 0
            cli::cli_alert_info(
                stringr::str_interp("Parsing genbank file [${gb_file}]"))
            private$conn = file(gb_file, "r")
            while(TRUE) {
                line = readLines(private$conn, n = 1)
                if ( length(line) == 0 ) {
                    break
                }
                line_count <- line_count + 1
                private$tag_handler(line, line_count)
            }
            close(private$conn)
            cli::cli_alert_success(
                stringr::str_interp("Done! [${line_count}] lines parsed"))
        },

        in_tag = function() {
            return(!is.na(private$current_tag))
        },

        tag_handler = function(line, pos) {
            # if leading whitespace continue quickly
            if (grepl("^ +", line) && private$in_tag()) {
                if (
                    grepl("^     [^ ]", line) &&
                    private$current_tag == "FEATURES") {
                    private$flush_feature()
                    private$extend_tag(line)
                } else {

                    if (private$in_tag()) {
                        private$extend_tag(line)
                    }
                }
            } else {

                tagged <- FALSE
                for (tag in private$key_tags) {
                    if (grepl(paste0("^", tag), line)) {
                        private$enter_tag(tag, pos)
                        private$extend_tag(line)
                        tagged <- TRUE
                    }
                }
                if (!tagged) {
                    private$exit_tag(pos)
                }
            }
        },


        enter_tag = function(tag, pos) {
            if (!is.na(private$current_tag)) {
                private$exit_tag(pos)
            }
            # cli::cli_alert(
            #     stringr::str_interp("entering tag [${tag}] at line [${pos}]"))
            private$current_tag <- tag
        },


        exit_tag = function(pos) {
            if (!is.na(private$current_tag)) {
                # cli::cli_alert(
                #     stringr::str_interp(
                #         "leaving tag [${private$current_tag}] at line [${pos}]"))

                if (private$current_tag == "ACCESSION") {
                    self$accession <- private$clip_tag()
                } else if (private$current_tag == "VERSION") {
                    self$version <- private$clip_tag()
                } else if (private$current_tag == "DEFINITION") {
                    self$definition <- private$clip_tag()
                } else if (private$current_tag == "FEATURES") {
                    # features are messily parsed during load - there should
                    # be a single residual feature to process ...
                    private$flush_feature()
                } else if (private$current_tag == "ORIGIN") {
                    private$process_sequence()
                }

                private$current_tag <- NA
                private$tag_string <- NULL
            }
        },


        clip_tag = function() {
            regex = paste0("(?<=",private$current_tag,").+")
            return(
                stringr::str_trim(
                    stringr::str_extract(
                        private$tag_string, regex))
           )
        },


        extend_tag = function(line) {
            private$tag_string <- append(private$tag_string, line)
        },


        flush_feature = function(silent=TRUE) {
            #cli::cli_alert_warning("feature flush called")
            private$debug_counter <- private$debug_counter + 1
            if (!silent) {
                print(private$tag_string)
            }

            if (grepl("^     source", private$tag_string[1])) {
                private$feature_source()
            } else if (grepl("^     CDS", private$tag_string[1])) {
                private$feature_cds()
            }

            private$tag_string <- NULL
            if (private$debug_counter > 5) {
                # silent_stop("debugging")
            }
        },


        feature_source = function() {
            #print(private$tag_string)

            #silent_stop("debugging")
        },


        feature_cds = function(subfeature_tag = "     CDS +") {

            # clip the subfeature tag & leading whitespace
            private$tag_string <- private$tag_string %>%
                stringr::str_replace(subfeature_tag, "") %>%
                stringr::str_trim() %>%
                stringr::str_c(collapse="")
            # extract start nucleotide, end nucleotide and strand
            coordinates <- private$extract_positions(
                private$tag_string %>%
                    stringr::str_extract("^[^/]+"))

            gene_id <- private$tag_string %>%
                private$extract_value_from_regex(
                    "(?<=gene=\")[^\"]+", "(?<=locus_tag=\")[^\"]+")

            translation <- stringr::str_extract(
                private$tag_string, "(?<=translation=\")[^\"]+")

            #cli::cli_alert(
            #    stringr::str_interp("[${gene_id}] --> [${coordinates}]"))
            private$gb_cds <- private$gb_cds %>%
                dplyr::bind_rows(
                    list(
                        chr=coordinates[1],
                        gene=gene_id,
                        strand=coordinates[2],
                        start=as.integer(coordinates[3]),
                        end=as.integer(coordinates[4]),
                        translation=translation
                        )
                    )
        },


        extract_value_from_regex = function(text, regex_a, regex_b) {
            if (stringr::str_detect(text, regex_a)) {
                stringr::str_extract(text, regex_a)
            } else if (stringr::str_detect(text, regex_b)) {
                stringr::str_extract(text, regex_b)
            } else {
                silent_stop("broken")
            }
        },


        extract_positions = function(position_str) {
            strand <- "+"
            if (stringr::str_detect(position_str, "complement")) {
                strand <- "-"
                position_str <- stringr::str_extract(
                    position_str, "(?<=complement\\()[^\\)]+")
            }
            positions <- unlist(stringr::str_split(position_str, "\\.\\."))

            # there are a few edge cases where gene locations are not atomic
            filter_position = function(str) {
                str_mod <- str %>%
                    stringr::str_replace("<", "") %>%
                    stringr::str_replace(">", "")

                xval = stringr::str_extract(str_mod, "[^\\d]+")
                if (is.na(xval)) return(str_mod)

                print(str)
                silent_stop("numeric issue")
            }

            return(c(self$accession,
                     strand,
                     filter_position(positions[1]),
                     filter_position(positions[2])))
        },


        process_sequence = function() {
            cli::cli_alert(stringr::str_interp("processing sequence [${length(private$tag_string)}] lines"))
        }

    )


)
