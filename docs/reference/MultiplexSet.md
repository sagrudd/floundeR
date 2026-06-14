# R6 Class for loading, visualising and analysing barcode information

R6 Class for loading, visualising and analysing barcode information

R6 Class for loading, visualising and analysing barcode information

## Super class

[`floundeR::FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
-\> `MultiplexSet`

## Active bindings

- `enumerate`:

  prepares a simple `2D Angenieux` enumeration of the provided dataset
  for quick visualisation of the dataset.

## Methods

### Public methods

- [`MultiplexSet$new()`](#method-MultiplexSet-new)

- [`MultiplexSet$as_tibble()`](#method-MultiplexSet-as_tibble)

- [`MultiplexSet$read_length_bins()`](#method-MultiplexSet-read_length_bins)

- [`MultiplexSet$quality_bins()`](#method-MultiplexSet-quality_bins)

- [`MultiplexSet$sequencingset()`](#method-MultiplexSet-sequencingset)

- [`MultiplexSet$temporalset()`](#method-MultiplexSet-temporalset)

- [`MultiplexSet$clone()`](#method-MultiplexSet-clone)

Inherited methods

- [`floundeR::FloundeR$bin_data()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-bin_data)
- [`floundeR::FloundeR$num_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-num_scale)
- [`floundeR::FloundeR$nums_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-nums_scale)
- [`floundeR::FloundeR$print()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Initialise a new instance of the R6 Class `MultiplexSet`

#### Usage

    MultiplexSet$new(seqsum = NA, barcoding_summary_file = NA)

#### Arguments

- `seqsum`:

  a tibble of sequencing summary information

- `barcoding_summary_file`:

  is a file.path to the corresponding barcoding_summary file that should
  be merged with the `seqsum` content.

------------------------------------------------------------------------

### Method `as_tibble()`

Export the imported dataset(s) as a tibble

This object consumes a sequencing summary file (and optionally the
corresponding barcoding_summary file) and creates an object in memory
that can be explored, sliced and filtered. This method dumps out the
in-memory object for further exploration and development.

#### Usage

    MultiplexSet$as_tibble()

#### Returns

A tibble representation of the starting dataset

------------------------------------------------------------------------

### Method `read_length_bins()`

bin the sequences in `seqsum` content into bins of sequence length

The nanopore sequencing run is expected to return a collection of
sequences that vary in their length distributions; this variance is a
function of the sequencing library prepared, the starting DNA etc. This
method is used to bin reads into uniform bins & assess the distribution
of sequence lengths.

#### Usage

    MultiplexSet$read_length_bins(
      qfilt = TRUE,
      normalised = TRUE,
      cumulative = FALSE,
      bins = 20,
      outliers = 0.025
    )

#### Arguments

- `qfilt`:

  - specifies how the quality information should be filtered

  - at the moment this only defines whether reporting is based on PASS
    or FAIL reads; would make more sense to have filtered by PASS / ALL?

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

    MultiplexSet$quality_bins(bins = 20, outliers = 0)

#### Arguments

- `bins`:

  the number of sequence bins that should be prepared (20 by default)

- `outliers`:

  defines the number of outliers (0 = 0%) that are excluded from the
  reads - should probably be deprecated for simplicity?

#### Returns

Angenieux 2D graph object

------------------------------------------------------------------------

### Method `sequencingset()`

Prepare a SequencingSet object from a given barcode

This method is used to subset the sequencing_summary information to
focus on a single barcode for more detailed analysis.

#### Usage

    MultiplexSet$sequencingset(barcode)

#### Arguments

- `barcode`:

  the barcode that should be reported

#### Returns

`SequencingSet` object containing information on barcode of interest.

------------------------------------------------------------------------

### Method `temporalset()`

Prepare a TemporalSet object from a given barcode

This method is used to subset the SequencingSummary information to focus
on a single barcode for a more detailed analysis of content.

#### Usage

    MultiplexSet$temporalset(barcode)

#### Arguments

- `barcode`:

  the barcode that should be reported

#### Returns

`TemporalSet` object containing information on barcode of interest.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MultiplexSet$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
