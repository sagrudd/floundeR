# R6 Class for loading and analysing nanopore sequencing_summary files

This class aims to simplify the handling and exploration of
sequencing_summary files and provides simple methods for accessing
information that can be used to assess the performance of a run.

## Super class

[`floundeR::FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
-\> `SequencingSummary`

## Public fields

- `sequencing_summary_file`:

  the file.path to the query sequencing summary file

- `barcoding_summary_file`:

  the file.path to the query barcoding summary file

## Active bindings

- `flowcell`:

  The sequencing summary file contains a collection of data that
  includes channel data. The channel data can be overlaid on spatial
  representations of flowcell layout to address spatial issues and to
  visualise the overall flowcell characteristics. This method creates a
  flowcell object that is suitable for these purposes.

- `sequencingset`:

  The `sequencingset` active binding returns a sequencingset object that
  is canonically structured around the `passes_filtering` logical field
  to allow assessment of sequencing characteristics.

- `temporalset`:

  The `temporalset` active binding prepares a temporalset object that is
  suitable for the temporal analysis of information within the
  sequencing summary file.

- `demultiplex`:

  The `demultiplex` active binding prepared a multiplexset object that
  can be used to explore the barcoded content contained within the
  sequencing summary file and to access attributes that are related to
  these information.

## Methods

### Public methods

- [`SequencingSummary$new()`](#method-SequencingSummary-new)

- [`SequencingSummary$as_tibble()`](#method-SequencingSummary-as_tibble)

- [`SequencingSummary$clone()`](#method-SequencingSummary-clone)

Inherited methods

- [`floundeR::FloundeR$bin_data()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-bin_data)
- [`floundeR::FloundeR$num_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-num_scale)
- [`floundeR::FloundeR$nums_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-nums_scale)
- [`floundeR::FloundeR$print()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Creates a new SequencingSummary object. This initialisation method
performs other sanity checking of the defined file(s) to ensure that it
is indeed parseable and creates the required data structures.

#### Usage

    SequencingSummary$new(sequencing_summary_file, barcoding_summary_file = NA)

#### Arguments

- `sequencing_summary_file`:

  The source sequencing_summary file.

- `barcoding_summary_file`:

  The source barcoding_summary_file file.

#### Returns

A new `SequencingSummary` object.

#### Examples

    sequencing_summary <- flnDr("sequencing_summary.txt.bz2")
    barcodes_summary <- flnDr("barcoding_summary.txt.bz2")
    seqsum <- SequencingSummary$new(sequencing_summary, barcodes_summary)

------------------------------------------------------------------------

### Method `as_tibble()`

Export the imported dataset(s) as a tibble

This object consumes a sequencing summary file (and optionally the
corresponding barcoding_summary file) and creates an object in memory
that can be explored, sliced and filtered. This method dumps out the
in-memory object for further exploration and development.

#### Usage

    SequencingSummary$as_tibble()

#### Returns

A tibble representation of the starting dataset

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SequencingSummary$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r

## ------------------------------------------------
## Method `SequencingSummary$new`
## ------------------------------------------------

sequencing_summary <- flnDr("sequencing_summary.txt.bz2")
barcodes_summary <- flnDr("barcoding_summary.txt.bz2")
seqsum <- SequencingSummary$new(sequencing_summary, barcodes_summary)
```
