
#' R6 Class for loading and analysing sequence sets
#'
#' @importFrom tidyr drop_na
#' @import emojifont
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
Infographic <- R6::R6Class(
    inherit = FloundeR,
    classname = "Infographic",
    public = list(
        #' @field panel.width defines width of an infographic block
        panel.width = 6,
        #' @field panel.height defines height of an infographic block
        panel.height = 4,
        #' @field panel.spacer defines spacing around an infographic block
        panel.spacer = 0.5,
        #' @field panel.x.offset defines x offset spacing for whole graphic
        panel.x.offset = 2,
        #' @field panel.y.offset defines y offset spacing for whole graphic
        panel.y.offset = 2,
        #' @field columns defines the number of columns to use in infographic
        #' plot
        columns = 4,


        #' @description
        #' Initialise a new instance of the R6 Class `Infographic`
        initialize = function() {
            private$.plot_elements <- list()
        },


        #' @description
        #' add an `InfographicItem` to the Infographic plot.
        #'
        #' @param item an `InfographicItem`
        add = function(item) {
            if (!class(item)[1] == "InfographicItem") {
                stop("Can only add [InfographicItem] elements")
            }
            private$.plot_elements <- append(private$.plot_elements, item)
        },

        #' @description
        #' Export the contained `Infographic` dataset(s) as a tibble
        #'
        #' @return A tibble representation for all the data
        as_tibble = function() {
            figures <- length(private$.plot_elements)
            figure_x <- seq(figures)
            suppressWarnings(
                length(figure_x) <- prod(
                    dim(matrix(figure_x, ncol = self$columns))))
            pmat <- matrix(figure_x, ncol = self$columns, byrow = TRUE)

            extracts_coords <- function(x) {
                where <- which(pmat==x, arr.ind=TRUE)
                x = self$panel.x.offset +
                    ((where[2]-1) * (self$panel.width + self$panel.spacer))
                y = self$panel.y.offset +
                    ((where[1]-1) * (self$panel.height + self$panel.spacer))
                y <- y * -1
                return(c(x=x, y=y, h=self$panel.height, w=self$panel.width))
            }

            df <- tibble::as_tibble(
                do.call(rbind,
                        lapply(seq(figures), extracts_coords)),
                .name_repair="universal")
            df$y <- df$y + (min(df$y) * -1 + self$panel.y.offset)

            df$key <- unlist(
                lapply(private$.plot_elements, function(x){return(x$.key)}))
            df$value <- unlist(
                lapply(private$.plot_elements, function(x){return(x$.value)}))
            df$icon <- unlist(
                lapply(private$.plot_elements, function(x){return(x$.icon)}))
            df$colour <- rep("steelblue", figures)
            return(df)
        },

        #' @description
        #' Plot the infographic to file (and display it immediately)
        #'
        #' @param display_file the file to write to the infographic to (a temp
        #' file will be created and used by default).
        plot = function(
            display_file = tempfile(
                pattern="file", tmpdir=tempdir(), fileext=".png")) {

            df <- self$as_tibble()

            plot <- ggplot(
                df,
                aes_string(
                    "x", "y", height="h", width="w",
                    label="key", fill="colour")) +
                geom_tile(fill = private$.tile_bg) +
                geom_text(
                    color = private$.txt_key_colour, hjust="left",
                    nudge_y=-1.5, nudge_x=-2.6, size=5) +
                geom_text(
                    label = df$icon, family = "fontawesome-webfont",
                    colour = private$.icon_colour, size = 23, hjust = "right",
                    nudge_x = 2.85,nudge_y = 0.8) +
                geom_text(
                    label = df$value, size = 10,
                    color = private$.txt_value_colour, fontface = "bold",
                    nudge_x = -2.6, hjust = "left")  +
                coord_fixed() +
                scale_fill_brewer(type = "qual", palette =  "Dark2") +
                theme_void() + guides(fill = FALSE)

            save_x = (max(df$x)+panel.width+panel.spacer) * 0.6
            save_y = (max(df$y)+panel.height+panel.spacer)  * 0.6

            ggplot2::ggsave(
                display_file, plot = plot, device = "png", units = "cm",
                width = save_x, height = save_y, dpi = 180)
            plot(magick::image_read(display_file))

        },


        #' @description
        #' Display a collection of fontawesome based infographics for picking.
        #'
        #' The `fontawesome` collection of icons contains over 700 icons of
        #' which some are more useful / desirable than others. This accessory
        #' method is used to render an Infographic report that summarised the
        #' available icons within a predefined range - the intention here is to
        #' make the selection of fonts to use in infographics a little simpler
        #' and easier. This replaces a dodgy notebook approach that was used
        #' previously.
        #'
        #' @param file - a file.path to use to write the infographic to
        #' @param offset - an integer offset defining where we should start
        #' rendering from in a broad sequence.
        #' @param rows - the number of rows to fill with sequential data.
        #' @param columns - the corresponding number of columns.
        display_fa = function(file, offset=0, rows=10, columns=6) {
            fonts <- emojifont::search_fontawesome("")
            if (offset >= length(fonts)) {
                stop(
                    paste0(
                        "There are only [",length(fonts),
                        "] fonts = nonsense"))
            }
            fmax = min(length(fonts), (rows * columns) + offset)
            print_fonts <- fonts[seq.int(from=(offset+1), to=fmax)]

            ig <- Infographic$new()
            ig$columns <- columns
            for (ff in print_fonts) {
                igi <- InfographicItem$new(ff, which(fonts==ff), ff)
                ig$add(igi)
            }
            ig$plot(file)

        }
    ),


    active = list(
        #' @field items
        #' return an integer describing the number of items that is contained
        #' within the `Infographic`.
        items = function() {
            return(length(private$.plot_elements))
        }
    ),

    private = list(
        .plot_elements = NULL,
        .tile_bg = RColorBrewer::brewer.pal(9, "Blues")[7],
        .txt_key_colour = RColorBrewer::brewer.pal(9, "Blues")[3],
        .icon_colour = RColorBrewer::brewer.pal(9, "Blues")[5],
        .txt_value_colour = RColorBrewer::brewer.pal(9, "Blues")[2]
    )
)





#' R6 Class for loading and analysing sequence sets
#'
#' @importFrom tidyr drop_na
#' @importFrom emojifont fontawesome
#'
#' @export
InfographicItem <- R6::R6Class(
    inherit = FloundeR,
    classname = "InfographicItem",
    public = list(
        #' @field .key the infographic key e.g. ReadCount
        .key = NULL,
        #' @field .value the element's value e.g. 42
        .value = NULL,
        #' @field .icon  the fa-awesome code to use for the cartoon display
        .icon = NULL,

        #' @description
        #' Initialise a new instance of the R6 Class `InfographicItem`
        #'
        #' This class is used to contain the information that is subsequently
        #' rendered by the `Infographic` class.
        #'
        #' @param key the infographic key e.g. ReadCount
        #' @param value the element's value e.g. 42
        #' @param icon the fa-awesome code to use for the cartoon display
        initialize = function(key=NA, value=NA, icon=NA) {
            fonts <- emojifont::search_fontawesome("")
            if (any(c(is.na(key), is.na(value), is.na(icon)))) {
                stop("InfographicItem requires key, value, icon")
            } else if (!icon %in% fonts) {
                stop(paste0("Specified FontAwesomeIcon [",icon,"] not found"))
            }
            self$.key <- key
            self$.value <- value
            self$.icon <- emojifont::fontawesome(icon)
        }
    )
)
