# R6 Class for visualising floundeR based datasets

This class aims to provide aplanat-like visualisation abstraction for
the floundeR framework

## Super class

[`floundeR::FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
-\> `Angenieux`

## Public fields

- `colourMap`:

  - default colourMap for plots requiring discrete colours

## Active bindings

- `data`:

  A method to dump out the stored data from an `Angenieux` object

## Methods

### Public methods

- [`Angenieux$new()`](#method-Angenieux-new)

- [`Angenieux$add()`](#method-Angenieux-add)

- [`Angenieux$plot()`](#method-Angenieux-plot)

- [`Angenieux$to_file()`](#method-Angenieux-to_file)

- [`Angenieux$set_title()`](#method-Angenieux-set_title)

- [`Angenieux$clone()`](#method-Angenieux-clone)

Inherited methods

- [`floundeR::FloundeR$bin_data()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-bin_data)
- [`floundeR::FloundeR$num_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-num_scale)
- [`floundeR::FloundeR$nums_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-nums_scale)
- [`floundeR::FloundeR$print()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Creates a new Angenieux object. This initialisation method performs
other sanity checking of the defined file(s) and creates the required
data structures.

#### Usage

    Angenieux$new(key, value)

#### Arguments

- `key`:

  the datatype that is being passed to the object

- `value`:

  the data that is being passed to the object

#### Returns

A new `Angenieux` object.

------------------------------------------------------------------------

### Method `add()`

add an `AngenieuxDecoration` to the plot.

#### Usage

    Angenieux$add(item)

#### Arguments

- `item`:

  an `AngenieuxDecoration`

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Prepare and present an Angenieux plot

#### Usage

    Angenieux$plot(...)

#### Arguments

- `...`:

  parameters passed on to downstream methods - please see examples for
  further examples as to how Angenieux plots can be customised using
  this approach.

------------------------------------------------------------------------

### Method `to_file()`

Specify that Angenieux plot should be saved to file

When working at the console an Angenieux plot may be plotted directly to
the console. When preparing reports through Rmarkdown or Pkgdown a more
logical saving of plots to a discrete file location may make more sense.
The method is used to instruct Angenieux that the plot should be saved
to a given location and with a given file format.

#### Usage

    Angenieux$to_file(
      target_file,
      width = 18,
      height = 12,
      units = "cm",
      dpi = "print"
    )

#### Arguments

- `target_file`:

  the file with extension e.g. `figure1.png`

- `width`:

  the width of figure to save (12 by default)

- `height`:

  the height of figure to save (7.5 by default)

- `units`:

  the unit to use for height and width (cm by default)

- `dpi`:

  the plot resolution (print/300 by default)

#### Returns

the original Angenieux object (self)

------------------------------------------------------------------------

### Method `set_title()`

Set the title used in the given Angenieux plot

#### Usage

    Angenieux$set_title(title)

#### Arguments

- `title`:

  - the title to use on the plot

#### Returns

the original Angenieux object (self)

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Angenieux$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
