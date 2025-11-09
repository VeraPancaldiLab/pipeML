# Merge Boruta Results

Merge results from multiple Boruta runs to identify robust feature
selections.

## Usage

``` r
merge_boruta_results(
  importance_values,
  decisions,
  file_name,
  iterations,
  threshold,
  return = TRUE
)
```

## Arguments

- importance_values:

  A list of data frames with feature importance values from each
  iteration.

- decisions:

  A list of character vectors with the decision labels from each
  iteration.

- file_name:

  A string used for naming the output plot file.

- iterations:

  Integer. The number of Boruta iterations performed.

- threshold:

  A numeric value between 0 and 1. Features labeled as 'Confirmed' or
  'Tentative' in more than `threshold * iterations` will be retained.

- return:

  Logical. Whether to save plots in the "Results/" directory.

## Value

A list containing:

- A vector of confirmed features.

- A vector of tentative features.

- A data frame with median importance values and final decisions.
