# R6 Class for loading and analysing sequence sets

R6 Class for loading and analysing sequence sets

R6 Class for loading and analysing sequence sets

## Super class

[`floundeR::FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
-\> `SequencingSet`

## Active bindings

- `enumerate`:

  prepares a simple `1D Angenieux` enumeration of the provided dataset
  for quick visualisation of the dataset.

- `N50`:

  Calculate and return the N50 value for passed quality sequence reads
  in the current `SequencingSet` object

- `mean`:

  Calculate and return the mean sequence length for passed quality reads
  in the `SequencingSet` object

## Methods

### Public methods

- [`SequencingSet$new()`](#method-SequencingSet-new)

- [`SequencingSet$as_tibble()`](#method-SequencingSet-as_tibble)

- [`SequencingSet$read_length_bins()`](#method-SequencingSet-read_length_bins)

- [`SequencingSet$quality_bins()`](#method-SequencingSet-quality_bins)

- [`SequencingSet$clone()`](#method-SequencingSet-clone)

Inherited methods

- [`floundeR::FloundeR$bin_data()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-bin_data)
- [`floundeR::FloundeR$num_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-num_scale)
- [`floundeR::FloundeR$nums_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-nums_scale)
- [`floundeR::FloundeR$print()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Initialise a new instance of the R6 Class `SequencingSet`

#### Usage

    SequencingSet$new(keycol, seqsum = NA)

#### Arguments

- `keycol`:

  a pointer to the column of interest in the `seqsum` to direct parsing
  and exploration of the file.

- `seqsum`:

  a tibble of sequencing summary information

------------------------------------------------------------------------

### Method `as_tibble()`

Export the imported dataset(s) as a tibble

This object consumes a sequencing summary file (and optionally the
corresponding barcoding_summary file) and creates an object in memory
that can be explored, sliced and filtered. This method dumps out the
in-memory object for further exploration and development.

#### Usage

    SequencingSet$as_tibble()

#### Returns

A tibble representation of the starting dataset

------------------------------------------------------------------------

### Method `read_length_bins()`

bin the sequences in `seqsum` content into bins of sequence length

The nanopore sequencing run is expected to return a collection of
sequences that vary in their length distributions; this variance is a
function of the sequencing library prepared, the starting DNA etc. This
method is used to bin reads into uniform bins to assess the distribution
of sequence lengths.

#### Usage

    SequencingSet$read_length_bins(
      normalised = TRUE,
      cumulative = FALSE,
      bins = 20,
      outliers = 0.025
    )

#### Arguments

- `normalised`:

  - should the sequence collection be reported to normalise for the
    number of sequence bases sequenced or the number of sequence reads -
    TRUE by default to normalise for sequenced bases.

- `cumulative`:

  defines whether cumulative sequence bases (reads) are reported per bin
  (FALSE by default).

- `bins`:

  the number of sequence bins that should be prepared (20 by default)

- `outliers`:

  defines the number of outliers (0.025 = 2.5%) that are excluded from
  the longest reads to prepare a richer distribution visulation - the
  plots can be bothered by the long tail of mini-whales.

#### Returns

Angenieux 2D graph object

------------------------------------------------------------------------

### Method `quality_bins()`

bin the sequences in `seqsum` content into bins of quality

The nanopore sequencing run is expected to return a collection of
sequences that vary in their quality distributions; this variance is a
function of the sequencing library prepared, the starting DNA etc. This
method is used to bin reads into uniform quality bins to assess the
overall quality of the run and to identify potential issues

#### Usage

    SequencingSet$quality_bins(bins = 20, outliers = 0)

#### Arguments

- `bins`:

  the number of sequence bins that should be prepared (20 by default)

- `outliers`:

  defines the number of outliers (0 = 0%) that are excluded from the
  reads - should probably be deprecated for simplicity??

#### Returns

Angenieux 2D graph object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SequencingSet$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
