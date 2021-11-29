
#' @importFrom ggplot2 %+replace%
angenieux_theme <- function() {

  ggplot2::theme_bw() %+replace%
    ggplot2::theme(

      # ggplot2::scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")),
      # ggplot2::scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
    )

}
