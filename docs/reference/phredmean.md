# calculate mean Phred scores from list of Q values

This method calculates the mean quality value from a list of q-values;
handles the challenge of relative log-scaling of the Q and challenges of
a more linear mean; this avoids an artificial inflation of phred-scores

## Usage

``` r
phredmean(q)
```

## Arguments

- q:

  is a vector of q values

## Value

mean phred scaled q-value

## Examples

``` r
mean(c(20,30,40))
#> [1] 30
phredmean(c(20,30,40))
#> [1] 24.31798
```
