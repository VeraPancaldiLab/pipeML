# data_example_survival

Example dataset for survival analysis. Uses the lung cancer dataset from
the `survival` package (complete cases only).

## Usage

``` r
data_example_survival
```

## Format

A data frame with 167 samples as rows and 10 columns: the survival time
in days (`time`), the event indicator (`status`: 1 = death, 0 =
censored) and 8 covariates.

## Source

[`survival::lung`](https://rdrr.io/pkg/survival/man/lung.html), with
`status` recoded from 1 = censored / 2 = dead to 0 = censored / 1 =
death.

## Examples

``` r
data(data_example_survival)
head(data_example_survival)
#>   inst time status age sex ph.ecog ph.karno pat.karno meal.cal wt.loss
#> 2    3  455      1  68   1       0       90        90     1225      15
#> 4    5  210      1  57   1       1       90        60     1150      11
#> 6   12 1022      0  74   1       1       50        80      513       0
#> 7    7  310      1  68   2       2       70        60      384      10
#> 8   11  361      1  71   2       2       60        80      538       1
#> 9    1  218      1  53   1       1       70        80      825      16
```
