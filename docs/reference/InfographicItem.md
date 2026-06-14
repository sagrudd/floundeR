# R6 Class for loading and analysing sequence sets

R6 Class for loading and analysing sequence sets

R6 Class for loading and analysing sequence sets

## Super class

[`floundeR::FloundeR`](https://sagrudd.github.io/floundeR/reference/FloundeR.md)
-\> `InfographicItem`

## Public fields

- `.key`:

  the infographic key e.g. ReadCount

- `.value`:

  the element's value e.g. 42

- `.icon`:

  the fa-awesome code to use for the cartoon display

## Methods

### Public methods

- [`InfographicItem$new()`](#method-InfographicItem-new)

- [`InfographicItem$clone()`](#method-InfographicItem-clone)

Inherited methods

- [`floundeR::FloundeR$bin_data()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-bin_data)
- [`floundeR::FloundeR$num_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-num_scale)
- [`floundeR::FloundeR$nums_scale()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-nums_scale)
- [`floundeR::FloundeR$print()`](https://sagrudd.github.io/floundeR/reference/FloundeR.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Initialise a new instance of the R6 Class `InfographicItem`

This class is used to contain the information that is subsequently
rendered by the `Infographic` class.

#### Usage

    InfographicItem$new(key = NA, value = NA, icon = NA)

#### Arguments

- `key`:

  the infographic key e.g. ReadCount

- `value`:

  the element's value e.g. 42

- `icon`:

  the fa-awesome code to use for the cartoon display

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    InfographicItem$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
