# R6 Class for loading and analysing FASTA files

This class aims to simplify the handling and exploration of FASTA files
and provides simple methods for accessing information that can be used
to assess bulk contents from a FASTA file - the analysis framework is
provided by Rsamtools (alone).

## Super class

[`floundeR::FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
-\> `Fasta`

## Active bindings

- `sequencingset`:

  The `sequencingset` active binding returns a sequencingset object that
  is canonically structured around the `passes_filtering` logical field
  to allow assessment of sequencing characteristics.

## Methods

### Public methods

- [`Fasta$new()`](#method-Fasta-new)

- [`Fasta$sequence_chunks()`](#method-Fasta-sequence_chunks)

- [`Fasta$get_sequence_chunk()`](#method-Fasta-get_sequence_chunk)

- [`Fasta$get_tibble_chunk()`](#method-Fasta-get_tibble_chunk)

- [`Fasta$get_index()`](#method-Fasta-get_index)

- [`Fasta$count()`](#method-Fasta-count)

- [`Fasta$as_tibble()`](#method-Fasta-as_tibble)

- [`Fasta$clone()`](#method-Fasta-clone)

Inherited methods

- [`floundeR::FloundeR$bin_data()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-bin_data)
- [`floundeR::FloundeR$num_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-num_scale)
- [`floundeR::FloundeR$nums_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-nums_scale)
- [`floundeR::FloundeR$print()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Creates a new Fasta object. This initialisation method performs other
sanity checking of the defined file(s) to ensure that it is indeed
parseable and creates the required data structures.

#### Usage

    Fasta$new(fasta_file)

#### Arguments

- `fasta_file`:

  The source sequencing_summary file.

#### Returns

A new `Fasta` object.

#### Examples

    canonical_fasta <- flnDr("cluster_cons.fasta.bgz")
    fasta <- Fasta$new(canonical_fasta)

------------------------------------------------------------------------

### Method `sequence_chunks()`

Split the fasta sequence file explored by the package into sequence
chunks for e.g. import into a relational database.

#### Usage

    Fasta$sequence_chunks(chunk_size = 1000)

#### Arguments

- `chunk_size`:

  The number of fasta entries that should be contained within a single
  chunk (default: 1000)

#### Returns

an invisible integer that defines the number of possible chunks; this
can for example be iterated over

------------------------------------------------------------------------

### Method `get_sequence_chunk()`

Get a chunk of fasta sequences from a larger monolithic file

#### Usage

    Fasta$get_sequence_chunk(id = 1)

#### Arguments

- `id`:

  the chunk (see `$sequence_chunks()`) to extract sequence for - this
  must be an integer that is \> 0 and \<= sequence_chunks.

#### Returns

DNAStringSet containing the fasta entries corresponding to the specified
sequence chunk.

------------------------------------------------------------------------

### Method `get_tibble_chunk()`

Get a chunk of fasta sequences from a larger monolithic file as a tibble

#### Usage

    Fasta$get_tibble_chunk(id)

#### Arguments

- `id`:

  the chunk (see `$sequence_chunks()`) to extract sequence for - this
  must be an integer that is \> 0 and \<= sequence_chunks.

#### Returns

`tibble` containing the fasta entries corresponding to the specified
sequence chunk.

------------------------------------------------------------------------

### Method `get_index()`

return the Rsamtools FASTA index

#### Usage

    Fasta$get_index()

#### Returns

GRanges object describing the fasta elements contained within the
sequence file

------------------------------------------------------------------------

### Method `count()`

return the number of sequence elements contained within the sequence
file specified

#### Usage

    Fasta$count()

#### Returns

integer of fasta entries in file

------------------------------------------------------------------------

### Method `as_tibble()`

Export the imported dataset(s) as a tibble

The Fasta R6 object consumes a fasta format file and creates an object
in memory that can be explored, sliced and filtered. This method dumps
out the in-memory object for further exploration and development.

#### Usage

    Fasta$as_tibble()

#### Returns

A tibble representation of the starting dataset

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Fasta$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r

## ------------------------------------------------
## Method `Fasta$new`
## ------------------------------------------------

canonical_fasta <- flnDr("cluster_cons.fasta.bgz")
fasta <- Fasta$new(canonical_fasta)
#> 
#> ── creating floundR::fasta with [cluster_cons.fasta.bgz] ───────────────────────
#> ℹ index for [cluster_cons.fasta.bgz] found
#> ℹ loading fasta index [cluster_cons.fasta.bgz.idx]
#> ✔ [18656] fasta entries parsed
```
