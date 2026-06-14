# R6 Class for loading and analysing nanopore (and other) FASTQ files

This class aims to simplify the handling and exploration of FASTQ files
and provides simple methods for accessing information that can be used
to assess the contents of a FASTQ file.

## Super class

[`floundeR::FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
-\> `Fastq`

## Public fields

- `fastq_file`:

  the file.path to the query FASTQ file

## Active bindings

- `sequencingset`:

  The `sequencingset` active binding returns a sequencingset object that
  is canonically structured around the `passes_filtering` logical field
  to allow assessment of sequencing characteristics.

## Methods

### Public methods

- [`Fastq$new()`](#method-Fastq-new)

- [`Fastq$as_tibble()`](#method-Fastq-as_tibble)

- [`Fastq$sequence_chunks()`](#method-Fastq-sequence_chunks)

- [`Fastq$get_sequence_chunk()`](#method-Fastq-get_sequence_chunk)

- [`Fastq$clone()`](#method-Fastq-clone)

Inherited methods

- [`floundeR::FloundeR$bin_data()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-bin_data)
- [`floundeR::FloundeR$num_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-num_scale)
- [`floundeR::FloundeR$nums_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-nums_scale)
- [`floundeR::FloundeR$print()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Creates a new Fastq object. This initialisation method performs other
sanity checking of the defined file(s) to ensure that it is indeed
parseable and creates the required data structures.

#### Usage

    Fastq$new(fastq_file)

#### Arguments

- `fastq_file`:

  The source sequencing_summary file.

#### Returns

A new `Fastq` object.

#### Examples

    canonical_fastq <- flnDr("example.fastq.gz")
    fastq <- Fastq$new(canonical_fastq)

------------------------------------------------------------------------

### Method `as_tibble()`

Export the imported dataset(s) as a tibble

This object consumes a sequencing summary file (and optionally the
corresponding barcoding_summary file) and creates an object in memory
that can be explored, sliced and filtered. This method dumps out the
in-memory object for further exploration and development.

#### Usage

    Fastq$as_tibble()

#### Returns

A tibble representation of the starting dataset

------------------------------------------------------------------------

### Method `sequence_chunks()`

Split the fastq sequence file explored by the package into sequence
chunks for e.g. import into a relational database.

#### Usage

    Fastq$sequence_chunks(chunk_size = 10000)

#### Arguments

- `chunk_size`:

  The number of fastq entries that should be contained within a single
  chunk (default: 10000)

#### Returns

an invisible integer that defines the number of possible chunks; this
can for example be iterated over

------------------------------------------------------------------------

### Method `get_sequence_chunk()`

Get a chunk of fastq sequences from a larger monolithic file. This
method can be called for up to `$sequence_chunks()` times or until NULL
results are returned.

#### Usage

    Fastq$get_sequence_chunk()

#### Returns

tibble containing the fastq entries corresponding to the available
sequence chunk.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Fastq$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r

## ------------------------------------------------
## Method `Fastq$new`
## ------------------------------------------------

canonical_fastq <- flnDr("example.fastq.gz")
fastq <- Fastq$new(canonical_fastq)
#> → opening fastq stream
```
