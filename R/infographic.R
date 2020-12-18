
#' R6 Class for loading and analysing sequence sets
#'
#' @description
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
        panel.width = 6,
        panel.height = 4,
        panel.spacer = 0.5,
        panel.x.offset = 2,
        panel.y.offset = 2,
        columns = 4,

        initialize = function() {
            private$.plot_elements <- list()
        },

        add = function(item) {
            if (!class(item)[1] == "InfographicItem") {
                stop("Can only add [InfographicItem] elements")
            }
            private$.plot_elements <- append(private$.plot_elements, item)
        },

        as_tibble = function() {
            figures <- length(private$.plot_elements)
            figure_x <- seq(figures)
            suppressWarnings(length(figure_x) <- prod(dim(matrix(figure_x, ncol = self$columns))))
            pmat <- matrix(figure_x, ncol = self$columns, byrow = TRUE)

            extracts_coords <- function(x) {
                where <- which(pmat==x, arr.ind=TRUE)
                x= self$panel.x.offset + ((where[2]-1) * (self$panel.width + self$panel.spacer))
                y= self$panel.y.offset + ((where[1]-1) * (self$panel.height + self$panel.spacer))
                y <- y * -1
                return(c(x=x, y=y, h=self$panel.height, w=self$panel.width))
            }

            df <- tibble::as_tibble(do.call(rbind, lapply(seq(figures), extracts_coords)), .name_repair="universal")
            df$y <- df$y + (min(df$y) * -1 + self$panel.y.offset)

            df$key <- unlist(lapply(private$.plot_elements, function(x){return(x$.key)}))
            df$value <- unlist(lapply(private$.plot_elements, function(x){return(x$.value)}))
            df$icon <- unlist(lapply(private$.plot_elements, function(x){return(x$.icon)}))
            df$colour <- rep("steelblue", figures)
            return(df)
        },

        plot = function() {

            df <- self$as_tibble()

            plot <- ggplot(
                df,
                aes_string("x", "y", height="h", width="w", label="key", fill="colour")) +
                geom_tile(fill = private$.tile_bg) +
                geom_text(color = private$.txt_key_colour, hjust="left", nudge_y=-1.5, nudge_x=-2.6, size=5) +
                geom_text(label = df$icon, family = "fontawesome-webfont", colour = private$.icon_colour, size = 23, hjust = "right", nudge_x = 2.85,nudge_y = 0.8) +
                geom_text(label = df$value, size = 10, color = private$.txt_value_colour, fontface = "bold", nudge_x = -2.6, hjust = "left")  +
                coord_fixed() +
                scale_fill_brewer(type = "qual", palette =  "Dark2") +
                theme_void() + guides(fill = FALSE)

            save_x = (max(df$x)+panel.width+panel.spacer) * 0.6
            save_y = (max(df$y)+panel.height+panel.spacer)  * 0.6

            display_file = tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")

            ggplot2::ggsave(display_file, plot = plot, device =
                                "png", units = "cm", width = save_x, height = save_y, dpi = 180)
            plot(magick::image_read(display_file))

        }
    ),

    active = list(
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
#' @description
#'
#' @importFrom tidyr drop_na
#' @importFrom emojifont fontawesome
#'
#' @export
InfographicItem <- R6::R6Class(
    inherit = FloundeR,
    classname = "InfographicItem",
    public = list(
        .key = NULL,
        .value = NULL,
        .icon = NULL,

        initialize = function(key=NA, value=NA, icon=NA) {
            if (any(c(is.na(key), is.na(value), is.na(icon)))) {
                stop("InfographicItem requires key, value, icon")
            }
            self$.key <- key
            self$.value <- value
            self$.icon <- emojifont::fontawesome(icon)
        }
    )
)
