% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/angenieux.R
\name{AngenieuxDecoration}
\alias{AngenieuxDecoration}
\title{R6 Class for describing additional Angenieux decorations}
\description{
R6 Class for describing additional Angenieux decorations

R6 Class for describing additional Angenieux decorations
}
\section{Super class}{
\code{\link[floundeR:FloundeR]{floundeR::FloundeR}} -> \code{AngenieuxDecoration}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{decoration}}{This field contains the decoration that will be applied to the Angenieux
object}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-AngenieuxDecoration-new}{\code{AngenieuxDecoration$new()}}
\item \href{#method-AngenieuxDecoration-.add_vline}{\code{AngenieuxDecoration$.add_vline()}}
\item \href{#method-AngenieuxDecoration-.add_vlegend}{\code{AngenieuxDecoration$.add_vlegend()}}
\item \href{#method-AngenieuxDecoration-.add_ggplot2}{\code{AngenieuxDecoration$.add_ggplot2()}}
\item \href{#method-AngenieuxDecoration-clone}{\code{AngenieuxDecoration$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="floundeR" data-topic="FloundeR" data-id="bin_data"><a href='../../floundeR/html/FloundeR.html#method-FloundeR-bin_data'><code>floundeR::FloundeR$bin_data()</code></a></span></li>
<li><span class="pkg-link" data-pkg="floundeR" data-topic="FloundeR" data-id="num_scale"><a href='../../floundeR/html/FloundeR.html#method-FloundeR-num_scale'><code>floundeR::FloundeR$num_scale()</code></a></span></li>
<li><span class="pkg-link" data-pkg="floundeR" data-topic="FloundeR" data-id="nums_scale"><a href='../../floundeR/html/FloundeR.html#method-FloundeR-nums_scale'><code>floundeR::FloundeR$nums_scale()</code></a></span></li>
<li><span class="pkg-link" data-pkg="floundeR" data-topic="FloundeR" data-id="print"><a href='../../floundeR/html/FloundeR.html#method-FloundeR-print'><code>floundeR::FloundeR$print()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-AngenieuxDecoration-new"></a>}}
\if{latex}{\out{\hypertarget{method-AngenieuxDecoration-new}{}}}
\subsection{Method \code{new()}}{
This is the constructor for the AngenieuxDecoration object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AngenieuxDecoration$new(decoration_type, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{decoration_type}}{This field is used to specify the type of decoration; the cleanest type
at the moment is currently the \code{ggplot2} type.}

\item{\code{...}}{the other variables passed on to methods contained within the
object}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-AngenieuxDecoration-.add_vline"></a>}}
\if{latex}{\out{\hypertarget{method-AngenieuxDecoration-.add_vline}{}}}
\subsection{Method \code{.add_vline()}}{
Add a vertical line to a ggplot2 graph within Angenieux
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AngenieuxDecoration$.add_vline(xintercept, colour = "green", size = 1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{xintercept}}{the point at which the vertical line will intercept the
x-axis}

\item{\code{colour}}{the colour of the line}

\item{\code{size}}{the width of the line (default 1)}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-AngenieuxDecoration-.add_vlegend"></a>}}
\if{latex}{\out{\hypertarget{method-AngenieuxDecoration-.add_vlegend}{}}}
\subsection{Method \code{.add_vlegend()}}{
Add a legend text to accompany a vertical line
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AngenieuxDecoration$.add_vlegend(
  xintercept,
  colour = "green",
  legend = "",
  hjust = 0,
  vjust = 1,
  size = 6
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{xintercept}}{the point at which the vertical line will intercept the
x-axis}

\item{\code{colour}}{the colour of the line}

\item{\code{legend}}{the text to display at the given location}

\item{\code{hjust}}{horizonal justify (0=left, 1=right)}

\item{\code{vjust}}{vertical justify (0=bottom, 1=top)}

\item{\code{size}}{the size of the font to present at the given location}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-AngenieuxDecoration-.add_ggplot2"></a>}}
\if{latex}{\out{\hypertarget{method-AngenieuxDecoration-.add_ggplot2}{}}}
\subsection{Method \code{.add_ggplot2()}}{
Just add some plain \code{ggplot2} to the AngenieuxDecoration and layer on to
the Angenieux plot - this is for the lazy hacking out and visualisation
of plots
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AngenieuxDecoration$.add_ggplot2(facet)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{facet}}{the stuff to be layered onto the graph}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-AngenieuxDecoration-clone"></a>}}
\if{latex}{\out{\hypertarget{method-AngenieuxDecoration-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AngenieuxDecoration$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
