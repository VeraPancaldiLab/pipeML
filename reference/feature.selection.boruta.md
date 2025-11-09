# Compute Feature Selection Using Repeated Boruta Algorithm

Repeatedly applies the Boruta feature selection algorithm and aggregates
results to determine consistently selected features.

## Usage

``` r
feature.selection.boruta(
  data,
  iterations = NULL,
  fix = FALSE,
  tentative = FALSE,
  doParallel = F,
  workers = NULL,
  file_name = NULL,
  threshold = NULL,
  return,
  verbose = T
)
```

## Arguments

- data:

  A data frame with the column "target" (factor) as the response and
  other columns as features.

- iterations:

  Integer. The number of Boruta iterations to perform.

- fix:

  Logical. If TRUE, applies TentativeRoughFix() to resolve tentative
  features after each iteration.

- tentative:

  Logical. Whether to include tentative features as confirmed in the
  training dataset.

- doParallel:

  Logical. Whether to use parallel processing.

- workers:

  Integer. Number of CPU cores to use for parallel execution. If NULL,
  uses all available cores minus one.

- file_name:

  A string for naming output plots and CSV files saved in the "Results/"
  directory.

- threshold:

  A numeric value between 0 and 1. A feature must be confirmed in more
  than `threshold * iterations` to be finally labeled as confirmed.

- return:

  Logical. Whether to save the resulting plots in the "Results/"
  directory.

- verbose:

  Boolen value to whether print or no the function messages

## Value

A list containing:

- A vector of confirmed features.

- A vector of tentative features.

- A data frame with median importance values and final decisions.
