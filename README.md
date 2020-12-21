
# floundeR

An R package to simplify the tidy analysis of Nanopore sequence data.

<!-- badges: start -->
[![Build Status](https://img.shields.io/github/languages/code-size/sagrudd/floundeR)](https://img.shields.io/github/languages/code-size/sagrudd/floundeR)
[![Releases](https://img.shields.io/github/downloads/sagrudd/floundeR/total)](https://img.shields.io/github/downloads/sagrudd/floundeR/total)
<!-- badges: end -->

The goal of floundeR is to provide a robust and Tidyverse compliant toolbox
that can be used to explore Nanopore sequence data. The analytical R code
contained within has been forked from earlier projects such as nanopoRe. The 
R code has been deduplicated and reimplemented as R6 objects. 

The project goals diverge from earlier versions of nanopoRe and the ambition is
now to enable abbreviated and simplified data exploration.

The documentation for the floundeR package can be found at
[https://sagrudd.github.io/floundeR/](https://sagrudd.github.io/floundeR/)

## Installation

There are a couple of dependencies for the successful usage of the `floundeR`
package.

* https://github.com/nanoporetech/vbz_compression/ needs to be installed on your
  computer if you are going to open FAST5 format files that have been compressed
  using the `VBZ` compression algorithm.

You can install the released version of floundeR from [github](
https://github.com/sagrudd) with:

``` r
install.packages("devtools")
devtools::install_github("sagrudd/floundeR")

# devtools does not appear to pull in some of the bioconductor stuff ...
install.packages("BiocManager")
BiocMananger::install("rhdf5")
# please see documents below - this requires the vbz plugin
```


## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(floundeR)
## basic example code
```

