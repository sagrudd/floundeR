# R6 Class for loading and analysing sequence sets

R6 Class for loading and analysing sequence sets

R6 Class for loading and analysing sequence sets

## Super class

[`floundeR::FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
-\> `Infographic`

## Public fields

- `panel.width`:

  defines width of an infographic block

- `panel.height`:

  defines height of an infographic block

- `panel.spacer`:

  defines spacing around an infographic block

- `panel.x.offset`:

  defines x offset spacing for whole graphic

- `panel.y.offset`:

  defines y offset spacing for whole graphic

- `columns`:

  defines the number of columns to use in infographic plot

## Active bindings

- `items`:

  return an integer describing the number of items that is contained
  within the `Infographic`.

## Methods

### Public methods

- [`Infographic$new()`](#method-Infographic-new)

- [`Infographic$add()`](#method-Infographic-add)

- [`Infographic$as_tibble()`](#method-Infographic-as_tibble)

- [`Infographic$plot()`](#method-Infographic-plot)

- [`Infographic$display_fa()`](#method-Infographic-display_fa)

- [`Infographic$clone()`](#method-Infographic-clone)

Inherited methods

- [`floundeR::FloundeR$bin_data()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-bin_data)
- [`floundeR::FloundeR$num_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-num_scale)
- [`floundeR::FloundeR$nums_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-nums_scale)
- [`floundeR::FloundeR$print()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Initialise a new instance of the R6 Class `Infographic`

#### Usage

    Infographic$new()

------------------------------------------------------------------------

### Method `add()`

add an `InfographicItem` to the Infographic plot.

#### Usage

    Infographic$add(item)

#### Arguments

- `item`:

  an `InfographicItem`

------------------------------------------------------------------------

### Method `as_tibble()`

Export the contained `Infographic` dataset(s) as a tibble

#### Usage

    Infographic$as_tibble()

#### Returns

A tibble representation for all the data

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Plot the infographic to file (and display it immediately)

#### Usage

    Infographic$plot(
      display_file = tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
    )

#### Arguments

- `display_file`:

  the file to write to the infographic to (a temp file will be created
  and used by default).

------------------------------------------------------------------------

### Method `display_fa()`

Display a collection of fontawesome based infographics for picking.

The `fontawesome` collection of icons contains over 700 icons of which
some are more useful / desirable than others. This accessory method is
used to render an Infographic report that summarised the available icons
within a predefined range - the intention here is to make the selection
of fonts to use in infographics a little simpler and easier. This
replaces a dodgy notebook approach that was used previously.

#### Usage

    Infographic$display_fa(file, offset = 0, rows = 10, columns = 6)

#### Arguments

- `file`:

  - a file.path to use to write the infographic to

- `offset`:

  - an integer offset defining where we should start rendering from in a
    broad sequence.

- `rows`:

  - the number of rows to fill with sequential data.

- `columns`:

  - the corresponding number of columns.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Infographic$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
