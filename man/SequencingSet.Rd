% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequencing_set.R
\name{SequencingSet}
\alias{SequencingSet}
\title{R6 Class for loading and analysing sequence sets}
\description{
R6 Class for loading and analysing sequence sets

R6 Class for loading and analysing sequence sets
}
\section{Super class}{
\code{\link[floundeR:FloundeR]{floundeR::FloundeR}} -> \code{SequencingSet}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{enumerate}}{prepares a simple \verb{1D Angenieux} enumeration of the provided dataset
for quick visualisation of the dataset.}

\item{\code{N50}}{Calculate and return the N50 value for passed quality sequence reads in
the current \code{SequencingSet} object}

\item{\code{mean}}{Calculate and return the mean sequence length for passed quality reads in
the \code{SequencingSet} object}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-SequencingSet-new}{\code{SequencingSet$new()}}
\item \href{#method-SequencingSet-as_tibble}{\code{SequencingSet$as_tibble()}}
\item \href{#method-SequencingSet-read_length_bins}{\code{SequencingSet$read_length_bins()}}
\item \href{#method-SequencingSet-quality_bins}{\code{SequencingSet$quality_bins()}}
\item \href{#method-SequencingSet-clone}{\code{SequencingSet$clone()}}
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
\if{html}{\out{<a id="method-SequencingSet-new"></a>}}
\if{latex}{\out{\hypertarget{method-SequencingSet-new}{}}}
\subsection{Method \code{new()}}{
Initialise a new instance of the R6 Class \code{SequencingSet}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SequencingSet$new(keycol, seqsum = NA)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{keycol}}{a pointer to the column of interest in the \code{seqsum} to
direct parsing and exploration of the file.}

\item{\code{seqsum}}{a tibble of sequencing summary information}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SequencingSet-as_tibble"></a>}}
\if{latex}{\out{\hypertarget{method-SequencingSet-as_tibble}{}}}
\subsection{Method \code{as_tibble()}}{
Export the imported dataset(s) as a tibble

This object consumes a sequencing summary file (and optionally the
corresponding barcoding_summary file) and creates an object in
memory that can be explored, sliced and filtered. This method dumps
out the in-memory object for further exploration and development.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SequencingSet$as_tibble()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A tibble representation of the starting dataset
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SequencingSet-read_length_bins"></a>}}
\if{latex}{\out{\hypertarget{method-SequencingSet-read_length_bins}{}}}
\subsection{Method \code{read_length_bins()}}{
bin the sequences in \code{seqsum} content into bins of sequence length

The nanopore sequencing run is expected to return a collection of
sequences that vary in their length distributions; this variance is a
function of the sequencing library prepared, the starting DNA etc. This
method is used to bin reads into uniform bins to assess the distribution
of sequence lengths.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SequencingSet$read_length_bins(
  normalised = TRUE,
  cumulative = FALSE,
  bins = 20,
  outliers = 0.025
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{normalised}}{\itemize{
\item should the sequence collection be reported to
normalise for the number of sequence bases sequenced or the number of
sequence reads - TRUE by default to normalise for sequenced bases.
}}

\item{\code{cumulative}}{defines whether cumulative sequence bases (reads) are
reported per bin (FALSE by default).}

\item{\code{bins}}{the number of sequence bins that should be prepared (20 by
default)}

\item{\code{outliers}}{defines the number of outliers (0.025 = 2.5\%) that are
excluded from the longest reads to prepare a richer distribution
visulation - the plots can be bothered by the long tail of mini-whales.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Angenieux 2D graph object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SequencingSet-quality_bins"></a>}}
\if{latex}{\out{\hypertarget{method-SequencingSet-quality_bins}{}}}
\subsection{Method \code{quality_bins()}}{
bin the sequences in \code{seqsum} content into bins of quality

The nanopore sequencing run is expected to return a collection of
sequences that vary in their quality distributions; this variance is a
function of the sequencing library prepared, the starting DNA etc. This
method is used to bin reads into uniform quality bins to assess the
overall quality of the run and to identify potential issues
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SequencingSet$quality_bins(bins = 20, outliers = 0)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{bins}}{the number of sequence bins that should be prepared (20 by
default)}

\item{\code{outliers}}{defines the number of outliers (0 = 0\%) that are
excluded from the reads - should probably be deprecated for simplicity??}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Angenieux 2D graph object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SequencingSet-clone"></a>}}
\if{latex}{\out{\hypertarget{method-SequencingSet-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SequencingSet$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
