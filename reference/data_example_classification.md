# data_example_classification

Example dataset for classification tasks. Derived from the Breast Cancer
Wisconsin dataset (`mlbench::BreastCancer`), cleaned by removing rows
with missing values.

## Usage

``` r
data_example_classification
```

## Format

A data frame with samples as rows and features as columns (e.g., cell or
tissue measurements).

## Source

`mlbench::BreastCancer` dataset

## Examples

``` r
data(data_example_classification)
head(data_example_classification)
#>   Cl.thickness Cell.size Cell.shape Marg.adhesion Epith.c.size Bare.nuclei
#> 1            5         1          1             1            2           1
#> 2            5         4          4             5            7          10
#> 3            3         1          1             1            2           2
#> 4            6         8          8             1            3           4
#> 5            4         1          1             3            2           1
#> 6            8        10         10             8            7          10
#>   Bl.cromatin Normal.nucleoli Mitoses target
#> 1           3               1       1      0
#> 2           3               2       1      0
#> 3           3               1       1      0
#> 4           3               7       1      0
#> 5           3               1       1      0
#> 6           9               7       1      1
```
