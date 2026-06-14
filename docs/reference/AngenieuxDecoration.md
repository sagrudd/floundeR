# R6 Class for describing additional Angenieux decorations

R6 Class for describing additional Angenieux decorations

R6 Class for describing additional Angenieux decorations

## Super class

[`floundeR::FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
-\> `AngenieuxDecoration`

## Public fields

- `decoration`:

  This field contains the decoration that will be applied to the
  Angenieux object

## Methods

### Public methods

- [`AngenieuxDecoration$new()`](#method-AngenieuxDecoration-new)

- [`AngenieuxDecoration$.add_vline()`](#method-AngenieuxDecoration-.add_vline)

- [`AngenieuxDecoration$.add_vlegend()`](#method-AngenieuxDecoration-.add_vlegend)

- [`AngenieuxDecoration$.add_ggplot2()`](#method-AngenieuxDecoration-.add_ggplot2)

- [`AngenieuxDecoration$clone()`](#method-AngenieuxDecoration-clone)

Inherited methods

- [`floundeR::FloundeR$bin_data()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-bin_data)
- [`floundeR::FloundeR$num_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-num_scale)
- [`floundeR::FloundeR$nums_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-nums_scale)
- [`floundeR::FloundeR$print()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-print)

------------------------------------------------------------------------

### Method `new()`

This is the constructor for the AngenieuxDecoration object

#### Usage

    AngenieuxDecoration$new(decoration_type, ...)

#### Arguments

- `decoration_type`:

  This field is used to specify the type of decoration; the cleanest
  type at the moment is currently the `ggplot2` type.

- `...`:

  the other variables passed on to methods contained within the object

------------------------------------------------------------------------

### Method `.add_vline()`

Add a vertical line to a ggplot2 graph within Angenieux

#### Usage

    AngenieuxDecoration$.add_vline(xintercept, colour = "green", size = 1)

#### Arguments

- `xintercept`:

  the point at which the vertical line will intercept the x-axis

- `colour`:

  the colour of the line

- `size`:

  the width of the line (default 1)

------------------------------------------------------------------------

### Method `.add_vlegend()`

Add a legend text to accompany a vertical line

#### Usage

    AngenieuxDecoration$.add_vlegend(
      xintercept,
      colour = "green",
      legend = "",
      hjust = 0,
      vjust = 1,
      size = 6
    )

#### Arguments

- `xintercept`:

  the point at which the vertical line will intercept the x-axis

- `colour`:

  the colour of the line

- `legend`:

  the text to display at the given location

- `hjust`:

  horizonal justify (0=left, 1=right)

- `vjust`:

  vertical justify (0=bottom, 1=top)

- `size`:

  the size of the font to present at the given location

------------------------------------------------------------------------

### Method `.add_ggplot2()`

Just add some plain `ggplot2` to the AngenieuxDecoration and layer on to
the Angenieux plot - this is for the lazy hacking out and visualisation
of plots

#### Usage

    AngenieuxDecoration$.add_ggplot2(facet)

#### Arguments

- `facet`:

  the stuff to be layered onto the graph

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    AngenieuxDecoration$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
