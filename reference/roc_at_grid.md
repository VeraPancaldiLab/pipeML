# Evaluate an ROC Curve on a Fixed False-Positive-Rate Grid (Internal)

Treats the ROC curve as a step function and returns, for each grid
value, the highest sensitivity reached at a false positive rate at or
below it.

## Usage

``` r
roc_at_grid(fpr, sensitivity, grid)
```

## Arguments

- fpr:

  Numeric vector of false positive rates, non-decreasing (as produced by
  [`get_sensitivity_specificity()`](https://verapancaldilab.github.io/pipeML/reference/get_sensitivity_specificity.md)).

- sensitivity:

  Numeric vector of sensitivities, same length as `fpr`.

- grid:

  Numeric vector of false positive rates in \[0, 1\].

## Value

Numeric vector of sensitivities, one per grid value (`NA` if the curve
is undefined, e.g. a resample without positives).
