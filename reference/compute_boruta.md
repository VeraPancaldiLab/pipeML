# Compute Boruta algorithm

Compute Boruta algorithm

## Usage

``` r
compute_boruta(data, seed, fix = TRUE)
```

## Arguments

- data:

  A data frame with the column of the variable to predict named "target"
  and the predictor features as additional columns.

- seed:

  A numeric value used to set the random seed for reproducibility.

- fix:

  Logical. If TRUE, applies TentativeRoughFix() from the Boruta package
  to resolve tentative features.

## Value

A list containing:

- A data frame with feature importance statistics.

- A character vector indicating the Boruta decision for each feature
  (Confirmed, Tentative, or Rejected).
