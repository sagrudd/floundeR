# Prepare a `flowcell` object from provided class

This method is intended to streamline the usage of floundeR in more
coherent workflows through the enabling of more magrittr-like piped
commands. This method will consume floundeR specified classes and return
a flowcell object where possible.

## Usage

``` r
to_flowcell(r6_object)
```

## Arguments

- r6_object:

  is the floundeR R6 object to extract flowcell from

## Value

Flounder flowcell object

## Examples

``` r
SequencingSummary$new(flnDr("sequencing_summary.txt.bz2")) %>% to_flowcell()
#> Preparing channel count information
#> <floundeR::Flowcell>
```
