
# floundeR

R package including accessory methods and toy datasets demonstrating the
analysis of Nanopore sequence data.

<!-- badges: start -->
[![Build Status](https://img.shields.io/github/languages/code-size/sagrudd/floundeR)](https://img.shields.io/github/languages/code-size/sagrudd/floundeR)
<!-- badges: end -->

The goal of floundeR is to provide a robust and Tidyverse compliant toolbox
that can be used to explore Nanopore sequence data. The analytical R code
contained within has been forked from earlier projects such as nanopoRe. The 
R code has been deduplicated and reimplemented as R6 objects. This diverges
from earlier versions of the code where there was a didactic emphasis - the
goal is now abbreviated and simplified data exploration.

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

