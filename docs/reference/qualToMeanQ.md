# calculate mean Phred score from an ASCII encoded phred string

FASTQ and BAM store per base qualities as an ASCII string. Accessory
methods in e.g. ShortRead allow for a sum of the numeric encoded scores;
this is not corrected for the log/linear so scores are synthetically
boosted

- this simple method performs mean on the character level data ...

## Usage

``` r
qualToMeanQ(qstr)
```

## Arguments

- qstr:

  is an ASCII encoded Phred quality score

## Value

mean phred scaled q-value

## Examples

``` r
qualToMeanQ('ABCDEF')
#> [1] 34.16954
```
