# R6 Class for analysing sequence sets with accompanying temporal data

R6 Class for analysing sequence sets with accompanying temporal data

R6 Class for analysing sequence sets with accompanying temporal data

## Super class

[`floundeR::FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
-\> `TemporalSet`

## Methods

### Public methods

- [`TemporalSet$new()`](#method-TemporalSet-new)

- [`TemporalSet$as_tibble()`](#method-TemporalSet-as_tibble)

- [`TemporalSet$runtime()`](#method-TemporalSet-runtime)

- [`TemporalSet$run_yield()`](#method-TemporalSet-run_yield)

- [`TemporalSet$feature_over_time()`](#method-TemporalSet-feature_over_time)

- [`TemporalSet$t50()`](#method-TemporalSet-t50)

- [`TemporalSet$clone()`](#method-TemporalSet-clone)

Inherited methods

- [`floundeR::FloundeR$bin_data()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-bin_data)
- [`floundeR::FloundeR$num_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-num_scale)
- [`floundeR::FloundeR$nums_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-nums_scale)
- [`floundeR::FloundeR$print()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Initialise a new instance of the R6 Class `TemporalSet`

#### Usage

    TemporalSet$new(seqsum = NA)

#### Arguments

- `seqsum`:

  a tibble of sequencing summary information

------------------------------------------------------------------------

### Method `as_tibble()`

Export the imported dataset(s) as a tibble

#### Usage

    TemporalSet$as_tibble()

#### Returns

A tibble representation of the starting dataset

------------------------------------------------------------------------

### Method `runtime()`

get a rounded runtime from temporal information in sequencing summary

The sequencing summary file contains information on pore dwell time and
start time for given sequencing reads. An analysis of the start time
information can be used to calculate a canonical sequencing run time
that is rounded to key temporal breaks.

#### Usage

    TemporalSet$runtime()

#### Returns

a numeric describing the rounded runtime in hours

------------------------------------------------------------------------

### Method `run_yield()`

Prepare summary statistics that describe a flowcell run's yield

A flowcell can yield both sequence bases and sequence reads and the
acquisition of these data has a temporal element. This method is used to
summarise run performance through assessment of yield per unit time.

#### Usage

    TemporalSet$run_yield(resolution = 15, bases = TRUE, cumulative = TRUE)

#### Arguments

- `resolution`:

  describes the temporal resolution (in minutes) by which yield will be
  summarised.

- `bases`:

  a logical that describes whether the method will summarise sequence
  reads or sequence bases - (TRUE) for bases by default

- `cumulative`:

  defines whether the number of reads(bases) per bin is described as
  actual number or as a temporally cumulative number

#### Returns

Angenieux object prepared for rendering in reports

------------------------------------------------------------------------

### Method `feature_over_time()`

Report a binned temporal data facet (such as speed, quality, length)

The temporal sequencing information can be used to summarise other data
facets in a time dependent manner to address questions such as whether
there is a change in sequencing speed, length or quality over time.

#### Usage

    TemporalSet$feature_over_time(
      resolution = 60,
      passes = TRUE,
      feature = "speed"
    )

#### Arguments

- `resolution`:

  describes the temporal resolution (in minutes) by which yield will be
  summarised (60 minutes by default).

- `passes`:

  is a logical that defines whether the plot should present data that
  has passed or failed QC - (TRUE) by default to select for only the QC
  passed sequence reads.

- `feature`:

  defines the feature to summarise (speed by default) but could include
  e.g. `sequence_length_template` or `mean_qscore_template`.

#### Returns

Angenieux object prepared for rendering in reports

------------------------------------------------------------------------

### Method `t50()`

calculate T50 timepoint at which 50% of sequence reads are obtained

Temporal data and yield data can be used to identify timepoints at which
a given fraction of the data has been obtained.

#### Usage

    TemporalSet$t50(passes = TRUE, t = 0.5)

#### Arguments

- `passes`:

  is a logical that defines that we should focus only on the passed QC
  reads (TRUE by default) - should probably be deprecated since asking
  for T50 of failed reads is just silly?

- `t`:

  is the fractional timepoint (0.5 by default) where we are interested
  in knowing the time at which this fraction of reads was obtained.

#### Returns

numeric describing timepoint in hours at which fractional data was
obtained.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    TemporalSet$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
