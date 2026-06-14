
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

You can install the released version of floundeR from [github](
https://github.com/sagrudd) with:

``` r
install.packages("devtools")
devtools::install_github("sagrudd/floundeR")
```

FAST5 support has been retired. The reboot direction for raw-signal QC is POD5
metadata and integrity inspection through Rust-backed package APIs.

Report rendering is moving to Grammateus semantic report contracts rather than
RMarkdown-first production reports. RMarkdown remains only for package
vignettes and transitional documentation; the boundary is documented in
[`LEGACY_REPORTING.md`](LEGACY_REPORTING.md).

Source installs now build a small embedded Rust extension. Install Cargo and
rustc first, then install the package from GitHub. Platform-specific source
install notes for macOS, Linux, Windows, Docker, and CI are recorded in
[`DEVELOPMENT.md`](DEVELOPMENT.md#source-install-requirements).

The detailed GitHub installation paths are recorded in
[`GITHUB_INSTALLATION.md`](GITHUB_INSTALLATION.md): public users install the
open-source package without private Grammateus assets, while authorized users
can configure optional prebuilt Grammateus runtimes for governed HTML/PDF
report rendering.

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
