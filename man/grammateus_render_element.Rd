\name{grammateus_render_element}
\alias{grammateus_render_element}
\alias{grammateus_render_figure_html}
\alias{grammateus_render_figure_pdf}
\title{Render Grammateus report elements through the Rust binding}
\usage{
grammateus_render_element(element, format = c("html", "pdf"))

grammateus_render_figure_html(element)

grammateus_render_figure_pdf(element)
}
\arguments{
\item{element}{A Grammateus semantic element prepared by floundeR, currently
a \code{flounder_grammateus_figure}.}

\item{format}{Output format: \code{html} or \code{pdf}.}
}
\value{
Rendered HTML as a character scalar for \code{format = "html"}, or PDF
bytes as a raw vector for \code{format = "pdf"}.
}
\description{
\code{grammateus_render_element()} sends a prepared floundeR Grammateus
semantic element to the compiled Rust report-rendering boundary.
\code{format = "html"} returns rendered HTML as a character scalar, while
\code{format = "pdf"} returns PDF bytes as a raw vector. Public floundeR builds
do not bundle the private Grammateus runtime; when it is not linked these
functions fail with a typed
\code{flounder_grammateus_runtime_unavailable} condition.
}
