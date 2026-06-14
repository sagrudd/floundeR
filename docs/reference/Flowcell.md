# R6 Class for performing Flowcell centric analyses

This class aims to simplify the handling and exploration of Flowcell
based data and contains various presets, designs and visualisation tools
required for assessing flowcell performance and metrics.

## Super class

[`floundeR::FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
-\> `Flowcell`

## Active bindings

- `platform`:

  Have a guess at the most likely flowcell platform used

  The sequencing summary file contains no information on the sequencing
  device or flowcell used. For the preparation of channel density maps
  it is worth considering which flowcell type is most likely to have
  been used - this can be guessed on the number of channels described
  within the data

- `density_data`:

  produce channelMap for spatial plots

  prepares a matrix of X, Y coordinates and the corresponding readcount
  information for the type of flowcell predicted by
  `get_flowcell_platform`

## Methods

### Public methods

- [`Flowcell$new()`](#method-Flowcell-new)

- [`Flowcell$set_channel_counts()`](#method-Flowcell-set_channel_counts)

- [`Flowcell$as_tibble()`](#method-Flowcell-as_tibble)

- [`Flowcell$clone()`](#method-Flowcell-clone)

Inherited methods

- [`floundeR::FloundeR$bin_data()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-bin_data)
- [`floundeR::FloundeR$num_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-num_scale)
- [`floundeR::FloundeR$nums_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-nums_scale)
- [`floundeR::FloundeR$print()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Creates a new Flowcell object. This initialisation method performs other
sanity checking of the defined file(s) and creates the required data
structures.

#### Usage

    Flowcell$new()

#### Returns

A new `Flowcell` object.

------------------------------------------------------------------------

### Method `set_channel_counts()`

set channel count summary information

This method is used to provide primitive channel count information for
the number of total reads that have been observed per channel - this is
used for the generation of spatial plots

#### Usage

    Flowcell$set_channel_counts(channel_counts)

#### Arguments

- `channel_counts`:

  a tibble of count information

------------------------------------------------------------------------

### Method `as_tibble()`

Export the imported dataset(s) as a tibble

This object consumes a sequencing summary file (and optionally the
corresponding barcoding_summary file) and creates an object in memory
that can be explored, sliced and filtered. This method dumps out the
in-memory object for further exploration and development.

#### Usage

    Flowcell$as_tibble()

#### Returns

A tibble representation of the starting dataset

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Flowcell$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
