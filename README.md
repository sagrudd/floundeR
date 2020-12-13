
# floundeR

<!-- badges: start -->
<!-- badges: end -->

The goal of floundeR is to ...

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

## FUBAR

* some issues with HDF5 parsing on macOS


## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(floundeR)
## basic example code
```

