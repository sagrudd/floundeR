# Pipe operator

Re-export of `magrittr`'s pipe operator for legacy floundeR workflows.

## Usage

``` r
lhs %>% rhs
```

## Arguments

- lhs:

  A value or expression to pass to the right-hand side.

- rhs:

  A function call or expression that receives `lhs`.

## Value

The result of evaluating `rhs` with `lhs` supplied by the pipe.
