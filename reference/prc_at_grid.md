# Evaluate a Precision-Recall Curve on a Fixed Recall Grid (Internal)

For each grid value, returns the precision at the first threshold whose
recall reaches that value.

## Usage

``` r
prc_at_grid(recall, precision, grid)
```

## Arguments

- recall:

  Numeric vector of recall values, non-decreasing (as produced by
  [`get_sensitivity_specificity()`](https://verapancaldilab.github.io/pipeML/reference/get_sensitivity_specificity.md)).

- precision:

  Numeric vector of precision values, same length as `recall`.

- grid:

  Numeric vector of recall values in \[0, 1\].

## Value

Numeric vector of precisions, one per grid value (`NA` if the curve is
undefined, e.g. a resample without positives).
