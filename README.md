
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


## `floundeR` and a BasicQC analysis to assess a flowcell run

BasicQC was the name of the original Nanopore tutorial that introduced an R
workflow to question a flowcell's performance on the basis of the 
`sequencing_summary` file produced by MinKNOW or Guppy (Albacore in the past).
With the transition to (Python based) [EPI2ME Labs](https://labs.epi2me.io) and
the deprecation of the original tutorial this workflow is being maintained to
support the requirements of R-aficionados.

Instead of using the packaged and toy dataset, let's use a slightly more
robust dataset to show what the tool can really do.

```
aws.s3::save_object(
   "/gm24385_2020.11/flowcells/20201026_1645_6B_PAG07165_d42912aa/sequencing_summary_PAG07165_2dfda515.txt", 
   bucket = "s3://ont-open-data/", 
   region="eu-west-1")
```


``` r
library(floundeR)
sequencing_summary <- "sequencing_summary_PAG07165_2dfda515.txt"
seqsum <- SequencingSummary$new(sequencing_summary)
seqsum$flowcell$density_data$plot
```

<img src="https://raw.githubusercontent.com/sagrudd/floundeR/main/docs/articles/figure_1.png" />

