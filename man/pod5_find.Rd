% Please edit documentation in R/pod5.R
\name{pod5_find}
\alias{pod5_find}
\title{Discover local POD5-containing folders}
\usage{
pod5_find(path)
}
\arguments{
\item{path}{Character scalar. Directory to search.}
}
\value{
A data frame with one row per discovered directory and the columns
\code{path}, \code{pod5_file_count}, \code{total_bytes},
\code{oldest_modified_utc}, and \code{newest_modified_utc}.
}
\description{
\code{pod5_find()} recursively searches a local directory tree for folders that
directly contain one or more \code{.pod5} files. Discovery is implemented
through the in-process Rust extension and the curated \code{pod5-tools} library
API; it does not invoke an external command-line process.
}
