# Evaluate an ROC Curve on a Fixed False-Positive-Rate Grid (Internal)

Returns, for each grid value, the sensitivity of the ROC curve at that
false positive rate, with straight lines between consecutive points: the
curve as
[`get_curves()`](https://verapancaldilab.github.io/pipeML/reference/get_curves.md)
draws it
([`geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html))
and as
[`calculate_auroc()`](https://verapancaldilab.github.io/pipeML/reference/calculate_auroc.md)
integrates it (trapezoids). At a vertical step (several points with the
same false positive rate), the highest sensitivity is used. Tied
probabilities of both classes give diagonal segments, which a step
function would place below the drawn curve.

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
