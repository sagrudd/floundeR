# R6 Class for floundeR based analyses

This is an overarching class which other flounder classes are expected
to inherit from. The `floundeR` object contains methods and variables
that are maintained for other objects - there is probably not a sensible
usecase for actually creating a `floundeR` instance on its own.

## Methods

### Public methods

- [`FloundeR$new()`](#method-FloundeR-new)

- [`FloundeR$print()`](#method-FloundeR-print)

- [`FloundeR$bin_data()`](#method-FloundeR-bin_data)

- [`FloundeR$num_scale()`](#method-FloundeR-num_scale)

- [`FloundeR$nums_scale()`](#method-FloundeR-nums_scale)

- [`FloundeR$clone()`](#method-FloundeR-clone)

------------------------------------------------------------------------

### Method `new()`

Creates a new FloundeR object.

#### Usage

    FloundeR$new()

#### Returns

A new `FloundeR` object.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

this print method overrides the standard print function as included with
R6 objects - this is to better define what is contained within an object
and to provide better floundeR abstraction

#### Usage

    FloundeR$print(...)

#### Arguments

- `...`:

  additional stuff passed on

#### Returns

nothing (at present) - output to stdout

------------------------------------------------------------------------

### Method `bin_data()`

general method for binning continuous data into sensible bins - used by
various methods within the broader floundeR framework.

#### Usage

    FloundeR$bin_data(data, bins = 20, outliers = 0.025, qmax = NA)

#### Arguments

- `data`:

  - the dataset (series) to be binned

- `bins`:

  - the number of bins to prepare

- `outliers`:

  - the fraction of longest outlying features that should be excluded to
    simplify the plot.

- `qmax`:

  - a manual upper limit to use

#### Returns

bin assignments for `data` elements

------------------------------------------------------------------------

### Method `num_scale()`

scale numerics in bases into kilo/mega/giga etc measurements for e.g.
plots and other graphic visualisations

#### Usage

    FloundeR$num_scale(val)

#### Arguments

- `val`:

  a numeric (raw)

#### Returns

character representation scaled accordingly

------------------------------------------------------------------------

### Method `nums_scale()`

Apply the `num_scale` method across a multi-member vector of numerics

#### Usage

    FloundeR$nums_scale(vals)

#### Arguments

- `vals`:

  the vector of numerics

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    FloundeR$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
