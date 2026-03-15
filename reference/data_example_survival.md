# data_example_survival

Example dataset for survival analysis. Uses the lung cancer dataset from
the `survival` package.

## Usage

``` r
data_example_survival
```

## Format

A data frame with samples as rows and variables as columns (e.g.,
survival time, status, covariates)

## Source

[`survival::lung`](https://rdrr.io/pkg/survival/man/lung.html)

## Examples

``` r
data(data_example_survival)
head(data_example_survival)
#>   inst time status age sex ph.ecog ph.karno pat.karno meal.cal wt.loss
#> 2    3  455      2  68   1       0       90        90     1225      15
#> 4    5  210      2  57   1       1       90        60     1150      11
#> 6   12 1022      1  74   1       1       50        80      513       0
#> 7    7  310      2  68   2       2       70        60      384      10
#> 8   11  361      2  71   2       2       60        80      538       1
#> 9    1  218      2  53   1       1       70        80      825      16
```
