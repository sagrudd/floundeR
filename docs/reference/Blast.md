# R6 Class for loading and analysing BLAST results in basic `Pairwise` format

BLAST results are a fundamental unit of basic comparative genomics. This
R6 object has been implemented for the systematic exploration of BLAST
results within the scope of the Lodestar project.

## Super class

[`floundeR::FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
-\> `Blast`

## Public fields

- `blast_file`:

  the file.path to the query BLAST results file

## Methods

### Public methods

- [`Blast$new()`](#method-Blast-new)

- [`Blast$count()`](#method-Blast-count)

- [`Blast$clone()`](#method-Blast-clone)

Inherited methods

- [`floundeR::FloundeR$bin_data()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-bin_data)
- [`floundeR::FloundeR$num_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-num_scale)
- [`floundeR::FloundeR$nums_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-nums_scale)
- [`floundeR::FloundeR$print()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Creates a new Blast object. This initialisation method performs other
sanity checking of the defined file(s) to ensure that it is indeed
parseable and creates the required data structures.

#### Usage

    Blast$new(blast_file)

#### Arguments

- `blast_file`:

  The source sequencing_summary file.

#### Returns

A new `Blast` object.

#### Examples

    blast_results <- flnDr("drosophila_uniref100.blastx.txt")
    blast <- Blast$new(blast_file=blast_results)

------------------------------------------------------------------------

### Method `count()`

Return the number of BLAST results that are contained within the BLAST
file provided.

#### Usage

    Blast$count()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Blast$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r

## ------------------------------------------------
## Method `Blast$new`
## ------------------------------------------------

blast_results <- flnDr("drosophila_uniref100.blastx.txt")
blast <- Blast$new(blast_file=blast_results)
#> → Parsing BLAST result file [drosophila_uniref100.blastx.txt]
```
