# Compute Concordance Index with Bootstrap Confidence Interval (Internal)

Computes the C-index (concordance index) for survival predictions and
estimates a 95% confidence interval via bootstrap resampling.

## Usage

``` r
compute_cindex_ci(
  data,
  time_col = "time",
  event_col = "event",
  pred_col = ".pred",
  n_boot = 1000,
  seed = 123
)
```

## Arguments

- data:

  Data frame containing survival outcomes and predictions.

- time_col:

  Character. Name of the survival time column (default = `"time"`).

- event_col:

  Character. Name of the event indicator column (default = `"event"`).

- pred_col:

  Character. Name of the prediction column (default = `".pred"`).

- n_boot:

  Integer. Number of bootstrap iterations for CI estimation (default =
  1000).

- seed:

  Integer. Random seed for reproducibility (default = 123).

## Value

A tibble with:

- c_index:

  Observed C-index on the input data.

- CI_lower:

  Lower bound of 95% bootstrap confidence interval.

- CI_upper:

  Upper bound of 95% bootstrap confidence interval.
