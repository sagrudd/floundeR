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


        get_accession = function() {

        },


        get_description = function() {

        },


        get_cds = function() {

        },


        get_sequence = function() {

        }



    ),

    private = list(
        conn = NA,
        current_tag = NA,
        tag_string = NULL,
        key_tags = c("DEFINITION", "ACCESSION", "FEATURES"),
        debug_counter = 0,
        cds = NA,

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
            cli::cli_alert(
                stringr::str_interp("entering tag [${tag}] at line [${pos}]"))
            private$current_tag <- tag
        },


        exit_tag = function(pos) {
            if (!is.na(private$current_tag)) {
                cli::cli_alert(
                    stringr::str_interp(
                        "leaving tag [${private$current_tag}] at line [${pos}]"))
                # print(private$tag_string)
                private$current_tag <- NA
                private$tag_string <- NULL
            }
        },


        extend_tag = function(line) {
            private$tag_string <- append(private$tag_string, line)
        },


        flush_feature = function() {
            #cli::cli_alert_warning("feature flush called")
            private$debug_counter <- private$debug_counter + 1
            #print(private$tag_string)

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

            # and coerce this collection of data into a GRanges
            cli::cli_alert(
                stringr::str_interp("[${gene_id}] --> [${coordinates}]"))

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
            return(position_str)
        }

    )


)
