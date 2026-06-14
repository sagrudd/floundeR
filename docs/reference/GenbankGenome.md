# R6 Class for loading and analysing Genbank whole genome files

This class aims to implement a couple of trivial methods for wrangling
whole microbial genome information from Genbank files

## Super class

[`floundeR::FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
-\> `GenbankGenome`

## Methods

### Public methods

- [`GenbankGenome$new()`](#method-GenbankGenome-new)

- [`GenbankGenome$list_cds()`](#method-GenbankGenome-list_cds)

- [`GenbankGenome$get_cds()`](#method-GenbankGenome-get_cds)

- [`GenbankGenome$focus_cds()`](#method-GenbankGenome-focus_cds)

- [`GenbankGenome$focus_range()`](#method-GenbankGenome-focus_range)

- [`GenbankGenome$as_tibble()`](#method-GenbankGenome-as_tibble)

- [`GenbankGenome$whole_genome_focus()`](#method-GenbankGenome-whole_genome_focus)

- [`GenbankGenome$focus_nucleotide()`](#method-GenbankGenome-focus_nucleotide)

- [`GenbankGenome$focus_codon()`](#method-GenbankGenome-focus_codon)

- [`GenbankGenome$get_nucleotide()`](#method-GenbankGenome-get_nucleotide)

- [`GenbankGenome$unfocus()`](#method-GenbankGenome-unfocus)

- [`GenbankGenome$is_focused()`](#method-GenbankGenome-is_focused)

- [`GenbankGenome$get_focus()`](#method-GenbankGenome-get_focus)

- [`GenbankGenome$clone()`](#method-GenbankGenome-clone)

Inherited methods

- [`floundeR::FloundeR$bin_data()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-bin_data)
- [`floundeR::FloundeR$num_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-num_scale)
- [`floundeR::FloundeR$nums_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-nums_scale)
- [`floundeR::FloundeR$print()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Creates a new GenbankGenome object. This initialisation method performs
minimal sanity checking of the defined file(s)

#### Usage

    GenbankGenome$new(gb_file)

#### Arguments

- `gb_file`:

  The source sequencing_summary file.

#### Returns

A new `GenbankGenome` object.

#### Examples

    TB_reference = flnDr("TB_H37Rv.gb.gz")
    tb <- GenbankGenome$new(TB_reference)

------------------------------------------------------------------------

### Method `list_cds()`

exports the annotated genomic features as a GenomicRanges object. This
can be used for the review of content in a genomics context.

#### Usage

    GenbankGenome$list_cds()

------------------------------------------------------------------------

### Method `get_cds()`

Get the CDS information for one or more annotated features from the
genome of interest.

#### Usage

    GenbankGenome$get_cds(feature_id = "fusA1")

------------------------------------------------------------------------

### Method `focus_cds()`

#### Usage

    GenbankGenome$focus_cds(feature_id = "fusA1")

------------------------------------------------------------------------

### Method `focus_range()`

#### Usage

    GenbankGenome$focus_range(chromosome = NULL, start = NULL, end = NULL)

------------------------------------------------------------------------

### Method `as_tibble()`

#### Usage

    GenbankGenome$as_tibble()

------------------------------------------------------------------------

### Method `whole_genome_focus()`

#### Usage

    GenbankGenome$whole_genome_focus()

------------------------------------------------------------------------

### Method `focus_nucleotide()`

#### Usage

    GenbankGenome$focus_nucleotide(position)

------------------------------------------------------------------------

### Method `focus_codon()`

#### Usage

    GenbankGenome$focus_codon(codon)

------------------------------------------------------------------------

### Method `get_nucleotide()`

#### Usage

    GenbankGenome$get_nucleotide(position, strand)

------------------------------------------------------------------------

### Method `unfocus()`

#### Usage

    GenbankGenome$unfocus()

------------------------------------------------------------------------

### Method `is_focused()`

#### Usage

    GenbankGenome$is_focused()

------------------------------------------------------------------------

### Method `get_focus()`

#### Usage

    GenbankGenome$get_focus()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    GenbankGenome$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r

## ------------------------------------------------
## Method `GenbankGenome$new`
## ------------------------------------------------

TB_reference = flnDr("TB_H37Rv.gb.gz")
tb <- GenbankGenome$new(TB_reference)
#> ℹ Parsing genbank file [/private/var/folders/g5/kp9259xj7778b2kzcmym7jd00000gn/T/RtmpfBxQW1/temp_libpath10eeb6d0106e7/floundeR/extdata/TB_H37Rv.gb.gz]
#> → processing sequence [73527] lines
#> → setting sequence
#> ✔ Done! [213179] lines parsed
```
